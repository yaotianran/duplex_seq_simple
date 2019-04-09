#!/usr/bin/env python3
# 0.5rc
import os, sys, random, time, argparse, collections, itertools, shutil, re, math, multiprocessing
import pysam, HTSeq, common

import os.path as path
import numpy as np

from multiprocessing import Process, Value, Lock, Manager
from Bio import SeqIO
from Bio import pairwise2 as Pairwise
from common import readable
from  common import cprint

SAMPLE_NAME = sys.argv[1]
IS_DUAL_END_TAG = True # Do the reads have tags on both 5'end and 3'end ?
TAG_LENGTH = 6  # How long is each tag on one end ?
MINIMUM_GROUP_MEMBER_AMOUNT = 8  # minimum segments amount in the group, should be an even integer
MAXIMUM_GROUP_MEMBER_AMOUNT = 20000

TAG_FALSE_MISMATCH_DETECTION = True
MAXIMUM_GROUP_POSITION_DEVIATION = 1
MINIMUM_INNER_TAG_SCORE = 11
GROUP_POSITION_DEVIATION_FACTOR = 1
MINIMUM_TAG_SCORE = 11  # total score of all 4 tags ( 2 for segment, 2 for its mate)

CORRECTION_VOTING_RATIO = 0.75 # a float greater than 0.5 or None for 'simple majority method'
#PREFER_INSERTION = True   # if a family have a same amount of insertion as non-insertion at a certain position, should we correct them as true insertion ?

GENOME_REF = '~/db/gencode.GRCh38.p12/GRCh38.primary_assembly.genome'
#GENOME_REF = '~/sda1/gencode.GRCh38.p12/GRCh38.primary_assembly.genome'
#CPU_INT = math.ceil( common.get_cpu_int() * 0.75 )
CPU_INT = 10
CORRECTRED_FASTQ = '{SAMPLE_NAME}.corrected.fastq'.format(SAMPLE_NAME= SAMPLE_NAME)
VERBOSE = False

BWA = 'bwa'
SAMTOOLS = 'samtools'
PICARD = 'java -jar ~/bin/picard_2_8_11.jar'

#MAXIMUM_POSITION_DEVIATION = 1   # maximum mapping position difference for all reads in a group/family
#MAXIMUM_LENGTH_DEVIATION = 1   # maximum read length difference for all reads in a group/family
#corrected_fastq_name_lst = []
#unsequenced_segments_lst = []  # to store the names of those segments whose query_sequence can't not be extracted.
#uncorrected_segment_dict = {} # key is a read name, value is tag string such as 'YA:i:1'

def __choose_group(group_tag, group_canidates):
    '''
    Accord segment and mate to choose which group_canidates keys they should belong to

    This function may raise an exception if an error occured

    Parameters:
        **group_tag**: namedtuple
            defined as collections.namedtuple('group_tag', ['segment_tag', 'mate_tag', 'chrom', 'pos']) # define the group tag

        **group_canidates**: list, in which every item is a group_tag
            its function

    Returns: a choosen_tag or None
        **choosen_tag**: a group_tag
            its meaning

        **None**
            doesn't belong to any group
    '''
    segment_tag_str = group_tag.segment_tag
    mate_tag_str = group_tag.mate_tag
    try:
        inner_score = Pairwise.align.globalms(segment_tag_str, mate_tag_str, 1,0,-1,-1 )[0][2]
        if inner_score < MINIMUM_INNER_TAG_SCORE:
            raise IndexError

    except IndexError:
        message = 'Inner score {inner_score}<{MINIMUM_INNER_TAG_SCORE}, too low. Segment tag: {segment_tag_str}. Mate tag: {mate_tag_str}.'.format(segment_tag_str= segment_tag_str,
                                                                                                                                                   mate_tag_str= mate_tag_str,
                                                                                                                                                   inner_score= inner_score,
                                                                                                                                                   MINIMUM_INNER_TAG_SCORE = MINIMUM_INNER_TAG_SCORE)
        if VERBOSE:
            print(message)
        return None

    total_tag_str = segment_tag_str + mate_tag_str
    map_position_float = group_tag.pos
    choosen_tag = None
    maximum_score = 0
    score_qualified_tags_lst = []
    for candidate_tag in group_canidates:
        if group_tag.chrom != candidate_tag.chrom or abs( map_position_float - candidate_tag.pos ) > MAXIMUM_GROUP_POSITION_DEVIATION:
            continue

        try:
            total_score = Pairwise.align.globalms(candidate_tag.segment_tag+candidate_tag.mate_tag, total_tag_str, 1,0,-1,-1 )[0][2] - abs(map_position_float - candidate_tag.pos)*GROUP_POSITION_DEVIATION_FACTOR
            if total_score < MINIMUM_TAG_SCORE:
                continue

            if total_score > maximum_score:
                maximum_score = total_score
                choosen_tag = candidate_tag
                break

        except IndexError:  # no hit found or anything make the current tag unqualified
            continue

    return choosen_tag


def __write_uncorrected_segment_to_fastq(SQL_dict, AlignedSegment, tag):
    '''
    Some segments are ungrouped or uncorrected for some reasons.

    Mark them by adding some special tags to their original fastq records and then write them to the end of UNCORRECTED_FASTQ file.

    Parameters:
    **SQL_dict**: Bio.File._SQLiteManySeqFilesDict object
        its function

    **AlignedSegment**: pysam.AlignedSegment object or string
        its function

    **tag**: string
        its function

    Returns:
        **1**: int
            {read_name} is not found in {SQL_dict}

        **2**: int
            Fail to write SeqRecord to UNCORRECTED_FASTQ file

        **3**: int
            Returncode of Bio.SeqIO.write is wrong.

        **4**: int
            Fail to get the correct query name
    '''
    # SQL_dict = SeqIO.index_db('temp/{}.idx'.format(SAMPLE_NAME), filenames=['temp/{}_1P_trimmed'.format(SAMPLE_NAME), 'temp/{}_2P_trimmed'.format(SAMPLE_NAME) ], format='fastq')

    # SeqRecord(seq=Seq('TTGGAGCTAGGTCCTTACTCTTCAGAAGGAGATAAAGGGGAAGGAAAGAATTTT...ATG', SingleLetterAlphabet()),
    #           id='E00500:287:HNVJLCCXY:3:1101:10003:24585/1',
    #           name='E00500:287:HNVJLCCXY:3:1101:10003:24585/1',
    #           description='E00500:287:HNVJLCCXY:3:1101:10003:24585/1 MI:Z:ATCG\tBC:Z:',
    #           dbxrefs=[])
    if isinstance(AlignedSegment, pysam.AlignedSegment):
        if AlignedSegment.is_read1:
            query_name = AlignedSegment.query_name + '/1'
        elif AlignedSegment.is_read2:
            query_name = AlignedSegment.query_name + '/2'
        else:
            return 4 # fail to get the correct query name
    elif isinstance(AlignedSegment, str):
        query_name = AlignedSegment

    else:
        return 4

    try:
        SeqRecord = SQL_dict[query_name]
    except:
        return 1

    SeqRecord.description += '\t' + tag

    try:
        with open(CORRECTRED_FASTQ, 'a') as ungroup_f:
            r = SeqIO.write(SeqRecord, ungroup_f, 'fastq')
            if r != 1:
                return 3
    except:
        return 2

    return 0
def __modify_tag(segment, tag_name: str) -> str:
    '''
    description='E00500:287:HNVJLCCXY:3:1101:10003:24585/1 Y5:Z:GGTTTT\tY3:Z:\tYL:i:151'
    __get_tag(description, 'Y5') will return 'GGTTTT'
    '''

    return

def __qualarray2qualstring(qual_array: np.array, offset= 33) -> bytes:
    '''
    Convert qual array to qual byte string
    qual array: array([38, 32, 38, 38, 38, 39, 39, 38, 38, 39, 39, 39, 38, 39, 38, 39, 39,
       37, 39, 39, 39, 38, 39, 39, 39, 39, 39, 38, 39, 37, 39, 39, 39, 39,
       39, 38, 37, 39, 39, 39, 39, 38, 32, 38, 39, 39, 37, 39, 39, 39, 39,
       33, 39, 38, 18, 39, 39, 39, 37, 35, 39, 38, 38, 38, 38, 38, 36, 32,
       36, 36, 38, 38, 38, 37, 36, 30, 30, 38, 36, 36, 19, 38, 38, 38, 37,
       37, 37, 37, 37, 37, 37, 33, 32, 32, 32, 32], dtype=uint8)

    qualstr: (b'GAACCAACCAAGCTCTCTTGAGGATCTTGAAGGAAACTGAATTCAAAAAGATCAAAGTGCTGGGCTCCGGTGCGTTCGGCACGGTGTATAAGGTAA')
    '''

    qualstr = b''
    for i in qual_array:
        qualstr += chr( int(i) + offset ).encode()

    return qualstr

def FastqPreprocess(fastq_file: str, is_forward: bool) -> None:
    '''
    Proprocess a fastq file, do the following progress:

    1. remove the 5'- and 3'- tags, write the to the tag

    2. add /1 or /2 to the end of read's name

    Parameters:
        **fastq_file**: type
            its function

        **is_forward**: type
            its function

    Returns:
        **value**: type
            its meaning
    '''
    #================format check================
    if not isinstance(fastq_file, str):
        message = 'fastq_file must be a string, currently is a {}'.format(type(fastq_file))
        raise TypeError(message)
    else:
        input_file_str = path.realpath(path.expanduser( fastq_file ))
        if not os.access(input_file_str, os.R_OK):
            message = '{} is not readable.'.format(fastq_file)
            raise FileNotFoundError(message)

    #================check over================
    output_file_str = 'temp/' + common.add_fix( path.split(input_file_str)[1], '_trimmed')
    output_f = open(output_file_str, 'w')

    tag_length = TAG_LENGTH
    i = 0
    for read in HTSeq.FastqReader(input_file_str):
        if (IS_DUAL_END_TAG and len(read) <= 2*tag_length) or (not IS_DUAL_END_TAG and len(read) <= tag_length):
            message = 'Length of read {name} (length={length}) is shorter than tag length {tag_length}. Skip'.format(name= read.name,
                                                                                                                     length= len(read),
                                                                                                                     tag_length= tag_length )
            if VERBOSE:
                common.cprint(message)
            continue

        # extract tag sequences
        try:
            orignal_name_str = read.name.split()[0]  # we only use the name before spaces. e.g. @M03639:59:000000000-G28YV:1:2102:14511:1891 1:N:0:TGTTGATT ----> @M03639:59:000000000-G28YV:1:2102:14511:1891
            if orignal_name_str[-2:] == '/1' or orignal_name_str[-2:] == '/2': # if the read's name end with '/1' or '/2'
                orignal_name_str = orignal_name_str[:-2]  # then remove it

            read_rc = read.get_reverse_complement()
            new_read = read
            if is_forward and IS_DUAL_END_TAG:  # if this is a forward-readthrough read
                tag_5_str = read.seq[:tag_length].decode()
                tag_3_str = read.seq[-tag_length:].decode()
                #new_read = read[tag_length:-tag_length]
                new_name = orignal_name_str + '/1'

            elif is_forward and not IS_DUAL_END_TAG:
                tag_5_str = read.seq[:tag_length].decode()
                tag_3_str = ''
                #new_read = read[tag_length:]
                new_name = orignal_name_str + '/1'

            elif not is_forward and IS_DUAL_END_TAG:
                tag_5_str = read_rc.seq[:tag_length].decode()
                tag_3_str = read_rc.seq[-tag_length:].decode()
                #new_read = read[tag_length:-tag_length]
                new_name = orignal_name_str + '/2'

            elif not is_forward and not IS_DUAL_END_TAG:
                tag_5_str = ''
                tag_3_str = read_rc.seq[:tag_length].decode()
                #new_read = read[tag_length:]
                new_name = orignal_name_str + '/2'

            else:
                pass

            new_read.name = '{new_name} Y5:Z:{tag_5_str}\tY3:Z:{tag_3_str}\tYL:i:{length}'.format(new_name= new_name,
                                                                                                  tag_3_str= tag_3_str,
                                                                                                  tag_5_str= tag_5_str,
                                                                                                  length= len(read) )
        except:
            message = 'Fail to extract tag sequences from read {}'.format(orignal_name_str)
            common.cprint(message)
            continue


        try:
            new_read.write_to_fastq_file(output_f)
            i += 1
            if i % 1000 == 0:
                print('{} sequences written.'.format(i), end='\r')
        except Exception as ex:
            common.cprint(ex)
            continue

    print('{} sequences written.'.format(i))
    output_f.close()

    return

def GroupSegments(samfile: str) -> dict:
    '''
    samfile must be sorted by querynames, and has unmapped and unpaired segments removed

    tag:
    YA:i:0
    YA:i:1 unmapped
    YA:i:2 supplementary_mapping
    YA:i:3 fail to get query sequence
    YA:i:4 reference_not_match
    YA:i:5 unpaired
    YA:i:6 improper_paired
    YA:i:7 corrupted_tags
    YA:i:8 mate_unmapped
    YA:i:9 too_few_members
    '''

    #===================format check===================
    if isinstance(samfile, str) and readable(samfile):
        samfile_str = path.realpath(path.expanduser(samfile))
    else:
        message = '{} is not a readable file.'.format(samfile)
        raise TypeError(message)
    #====================check over====================
    lock = Lock()
    group_tag = collections.namedtuple('group_tag', ['segment_tag', 'mate_tag', 'chrom', 'pos']) # define the group tag
    grouped_segments_dict = collections.defaultdict(list) # to store grouped segments {group_tag1:[segment1, segment2, segment3....], group_tag2:[segment4, segment5, segment6....]}
    grouped_segments_count_dict = collections.defaultdict(int) # to store how many segments in a group, {group_tag1:11, group_tag2:25, group_tag3:42.....}
    ungrouped_segments_count_dict = collections.defaultdict(int) # to store how many segments are not grouped
    grouped_segments_int = 0 # how many grouped segments
    ungrouped_segments_int = 0 # how many ungrouped segments

    original_read_dict = SeqIO.index_db('temp/{}.idx'.format(SAMPLE_NAME), filenames=['temp/{}_1P_trimmed'.format(SAMPLE_NAME), 'temp/{}_2P_trimmed'.format(SAMPLE_NAME) ], format='fastq')
    i = 0
    all_segments = pysam.AlignmentFile(samfile_str, 'r').fetch()
    for segment in all_segments:
        mate = next(all_segments)
        i += 1

        # deal with unpaired
        while segment.query_name != mate.query_name:
            ungrouped_segments_count_dict['unpaired'] += 1
            #r = __write_uncorrected_segment_to_fastq(original_read_dict, segment, 'YA:i:5')
            segment = mate
            mate = next(all_segments)

        # deal with wrong mapping position
        if segment.reference_name != mate.reference_name:
            ungrouped_segments_count_dict['reference_not_match'] += 1
            ungrouped_segments_int += 1
            r = __write_uncorrected_segment_to_fastq(original_read_dict, segment, 'YA:i:4')
            r = __write_uncorrected_segment_to_fastq(original_read_dict, mate, 'YA:i:4')
            continue


        # deal with improper_pair
        '''
        if not segment.is_proper_pair:
            ungrouped_segments_count_dict['improper_pair'] += 1
            r = __write_uncorrected_segment_to_fastq(original_read_dict, segment, 'YA:i:6')
            r = __write_uncorrected_segment_to_fastq(original_read_dict, mate, 'YA:i:6')
            continue
        '''

        #====================================================grouping======================================================
        segment_tag_str = segment.get_tag('Y5') + segment.get_tag('Y3')
        mate_tag_str = mate.get_tag('Y5') + mate.get_tag('Y3')
        pos = ( segment.reference_start + mate.reference_start ) / 2
        chrom = segment.reference_name
        current_tag = group_tag(segment_tag_str, mate_tag_str, chrom, pos)

        grouped_segments_dict[current_tag].extend([segment, mate])   # group the segment-pairs natively
        grouped_segments_count_dict[current_tag] += 1
        grouped_segments_int += 1

        if i % 10000 == 0:
            message = '[GROUP SEGMENTS] Processed segment pairs:{i}, Grouped pairs:{grouped_segments_int}, Ungrouped pairs:{ungrouped_segments_int} {stat}'.format(i= i,
                                                                                                                                                                grouped_segments_int= grouped_segments_int,
                                                                                                                                                                ungrouped_segments_int= ungrouped_segments_int,
                                                                                                                                                                stat= str(ungrouped_segments_count_dict)[27:-1])
            print(message, end='\r')

    message = '[GROUP SEGMENTS] Processed segment pairs:{i}, Grouped pairs:{grouped_segments_int}, Ungrouped pairs:{ungrouped_segments_int} {stat}'.format(i= i,
                                                                                                                                                        grouped_segments_int= grouped_segments_int,
                                                                                                                                                        ungrouped_segments_int= ungrouped_segments_int,
                                                                                                                                                        stat= str(ungrouped_segments_count_dict)[27:-1])
    print(message)

    if TAG_FALSE_MISMATCH_DETECTION:
        # pick up correct and wrong tags
        total = len(grouped_segments_dict)
        correct_tags_dict = collections.defaultdict(list)   # to store the correct group tags, {('chr1',19812413.0):[tag1, tag2], ('chr2',18365163.5):[tag3, tag4], .....}
        correct_tags_int = 0
        wrong_tags_dict = collections.defaultdict(list)
        wrong_tags_int = 0
        for tag in grouped_segments_dict.keys():
            if tag.segment_tag == tag.mate_tag and 'N' not in tag.segment_tag and 'N' not in tag.mate_tag: # this is a correct
                correct_tags_dict[ (tag.chrom, tag.pos) ].append(tag)
                correct_tags_int += 1
            else:
                wrong_tags_dict[ (tag.chrom, tag.pos) ].append(tag)
                wrong_tags_int += 1

        message = '[TAG_FALSE_MISMATCH_DETECTION] Total tag groups: {total}, correct groups: {correct_tags_int} at {i} map loci, wrong groups: {wrong_tags_int} at {j} map loci'.format(total           = total,
                                                                                                                                                                                    correct_tags_int= correct_tags_int,
                                                                                                                                                                                    i               = len(correct_tags_dict),
                                                                                                                                                                                    wrong_tags_int  = wrong_tags_int,
                                                                                                                                                                                    j               = len(wrong_tags_dict))
        print(message)

        # make the group correction
        t = time.time()
        i = 0
        for maploc_tu, wrong_tag_lst in wrong_tags_dict.items(): # to pick up qualified tag candidates according to map location. maploc_tu is a two-item tuple, such as ('chr1',19812413.0), wrong_tag_lst is a list containing a serial of group_tag
            candidate_tag_lst = []
            for mappos_float in np.arange(maploc_tu[1]- MAXIMUM_GROUP_POSITION_DEVIATION, maploc_tu[1]+MAXIMUM_GROUP_POSITION_DEVIATION+0.5, 0.5):
                candidate_tag_lst.extend( correct_tags_dict[ (maploc_tu[0], mappos_float) ] )
                if VERBOSE:
                    message = '[TAG_FALSE_MISMATCH_DETECTION] Added\t{candidates_int} candidate tags\tat {location}'.format(candidates_int= len( correct_tags_dict[(maploc_tu[0], mappos_float)] ),
                                                                              location= (maploc_tu[0], mappos_float) )
                    print(message)

            if VERBOSE:
                print('[TAG_FALSE_MISMATCH_DETECTION] Candidate tags: ', len(candidate_tag_lst) )

            for wrong_tag in wrong_tag_lst:
                if len(candidate_tag_lst) == 0:
                    similar_tag = None
                else:
                    similar_tag = __choose_group(wrong_tag, candidate_tag_lst)

                if similar_tag is None: # no similar group found
                    all_segments = iter(grouped_segments_dict[wrong_tag])
                    for segment in all_segments:
                        mate = next(all_segments)
                        ungrouped_segments_count_dict['corrupted_tags'] += 1
                        ungrouped_segments_int += 1
                        grouped_segments_int -= 1
                        r = __write_uncorrected_segment_to_fastq(original_read_dict, segment, 'YA:i:7')
                        r = __write_uncorrected_segment_to_fastq(original_read_dict, mate, 'YA:i:7')

                    grouped_segments_count_dict.pop(wrong_tag)
                    grouped_segments_dict.pop(wrong_tag)
                    message = '[TAG_FALSE_MISMATCH_DETECTION] No similar group found, Discard {wrong_tag}'.format(wrong_tag= wrong_tag)
                    if VERBOSE:
                        print(message)

                else:                   # we found a similar group and add the segment in the similar group and going to merge them together
                    member_amount_int = len(grouped_segments_dict[wrong_tag])
                    for segment in grouped_segments_dict.pop(wrong_tag):
                        if IS_DUAL_END_TAG and segment.is_read1:
                            segment.set_tag('Y5', similar_tag.segment_tag[:TAG_LENGTH], 'Z')
                            segment.set_tag('Y3', similar_tag.segment_tag[TAG_LENGTH:], 'Z')
                        elif IS_DUAL_END_TAG and segment.is_read2:
                            segment.set_tag('Y5', similar_tag.mate_tag[:TAG_LENGTH], 'Z')
                            segment.set_tag('Y3', similar_tag.mate_tag[TAG_LENGTH:], 'Z')
                        elif not IS_DUAL_END_TAG and segment.is_read1:
                            segment.set_tag('Y5', similar_tag.segment_tag[:TAG_LENGTH], 'Z')
                            segment.set_tag('Y3', similar_tag.segment_tag[TAG_LENGTH:], 'Z')
                        elif not IS_DUAL_END_TAG and segment.is_read2:
                            segment.set_tag('Y3', similar_tag.mate_tag[:TAG_LENGTH], 'Z')
                            segment.set_tag('Y5', similar_tag.mate_tag[TAG_LENGTH:], 'Z')
                        else:
                            raise TypeError
                        grouped_segments_dict[similar_tag].append(segment)

                    grouped_segments_count_dict[similar_tag] += grouped_segments_count_dict.pop(wrong_tag)
                    message = '[TAG_FALSE_MISMATCH_DETECTION] Merge {wrong_tag} (members: {member_amount_int}) into {similar_tag}.'.format(wrong_tag= wrong_tag,
                                                                                                                                           similar_tag= similar_tag,
                                                                                                                                           member_amount_int= member_amount_int)
                    if VERBOSE:
                        print(message)

                i += 1
                if i % 100 == 0:
                    message = '[TAG_FALSE_MISMATCH_DETECTION][Elapsed time: {time}s] Processed wrong groups: {i}/{wrong_tags_int}\t Ungrouped pairs: {ungroup} {stat}'.format(time= round(time.time()-t),
                                                                                                                                                                    i= i,
                                                                                                                                                                    wrong_tags_int= wrong_tags_int,
                                                                                                                                                                    ungroup= ungrouped_segments_int,
                                                                                                                                                                    stat= str(ungrouped_segments_count_dict)[27:-1])
                    print(message, end= '\r')

        message = '[TAG_FALSE_MISMATCH_DETECTION][Elapsed time: {time}s] Processed wrong groups: {i}/{wrong_tags_int}\t Ungrouped pairs: {ungroup} {stat}'.format(time= round(time.time()-t),
                                                                                                                                                        i= i,
                                                                                                                                                        wrong_tags_int= wrong_tags_int,
                                                                                                                                                        ungroup= ungrouped_segments_int,
                                                                                                                                                        stat= str(ungrouped_segments_count_dict)[27:-1])
        print(message)

        #====================================================group complete======================================================

    # discard the groups with too few members
    for tag, count in grouped_segments_count_dict.items():
        if count * 2 < MINIMUM_GROUP_MEMBER_AMOUNT:
            all_segments = iter(grouped_segments_dict[tag])
            for segment in all_segments:
                mate = next(all_segments)
                ungrouped_segments_count_dict['too_few_members'] += 1
                ungrouped_segments_int += 1
                grouped_segments_int -= 1
                r = __write_uncorrected_segment_to_fastq(original_read_dict, segment, 'YA:i:9')
                r = __write_uncorrected_segment_to_fastq(original_read_dict, mate, 'YA:i:9')

            #grouped_segments_count_dict.pop(tag)
            grouped_segments_dict.pop(tag)

            message = '[GROUP SEGMENTS] Dropping groups without enough member {MINIMUM_GROUP_MEMBER_AMOUNT}, Grouped pairs:{grouped_segments_int}, Ungrouped pairs:{ungrouped_segments_int} {stat}'.format(MINIMUM_GROUP_MEMBER_AMOUNT= MINIMUM_GROUP_MEMBER_AMOUNT,
                                                                                                                                                                                          grouped_segments_int= grouped_segments_int,
                                                                                                                                                                                          ungrouped_segments_int= ungrouped_segments_int,
                                                                                                                                                                                          stat= str(ungrouped_segments_count_dict)[27:-1])
            if grouped_segments_int % 1000 == 0:
                print(message, end= '\r')

    message = '[GROUP SEGMENTS] Dropping groups without enough member {MINIMUM_GROUP_MEMBER_AMOUNT}, Grouped pairs:{grouped_segments_int}, Ungrouped pairs:{ungrouped_segments_int} {stat}'.format(MINIMUM_GROUP_MEMBER_AMOUNT= MINIMUM_GROUP_MEMBER_AMOUNT,
                                                                                                                                                                                  grouped_segments_int= grouped_segments_int,
                                                                                                                                                                                  ungrouped_segments_int= ungrouped_segments_int,
                                                                                                                                                                                  stat= str(ungrouped_segments_count_dict)[27:-1])
    print(message)

    # do some statistics
    grouped_segments_count_lst = []
    for group_tag, pairs_lst in grouped_segments_dict.items():
        if len(pairs_lst) % 2 != 0:
            message = '{group_tag} contains uneven members {count}'.format(group_tag= group_tag, count= len(pairs_lst))
            cprint(message)
            continue
        grouped_segments_count_lst.append(len(pairs_lst))

    message = '[GROUP SEGMENTS] Member count statistics: total {total} tags\tMin. {min_int} Mean. {mean} Std. {std} Median {median} Max. {max_int}'.format(total= len(grouped_segments_dict),
                                                                                                                                           min_int= np.min(grouped_segments_count_lst),
                                                                                                                                           mean= round(np.mean(grouped_segments_count_lst),2),
                                                                                                                                           std= round(np.std(grouped_segments_count_lst),2),
                                                                                                                                           median= np.median(grouped_segments_count_lst),
                                                                                                                                           max_int= np.max(grouped_segments_count_lst))
    print(message)

    return grouped_segments_dict



def CorrectAndWrite(bamfile: str, output_file: str):
    '''
    Correct aligned bamfile and write every segment to a single interleaved fastq file. The name of each read will have a suffix of /1 or /2 to describe which end it came from.

    Banfile should be sorted by coordinate and indexed, in which all segments should belong to the same group and all paired

    Please note that he output fastq is *NOT* sorted (either by coordinate or by read name).

    Parameters:
        **bamfile**: string
            a bamfile, in which segments will be corrected, the segments should belong to the same group and all paired

        **lock**: multiprocessing.Lock object
            for locking the file handle if you use this function in multiprocessing

    Returns:
        **corrected_segments_amount**: int
            How many segments corrected in the group
    '''

    import multiprocessing
    from HTSeq import SequenceWithQualities as fastq

    #====================================================
    # to store all the SequenceWithQualities objects, so that durning the correction, we can correct the sequnces in the fastq_dict and not the sequnces of AlignedSegment objects themselves
    fastq_dict = collections.defaultdict(lambda : None)
    # to stote the [Y5, Y3, YL] tages for each segment
    tags_dict = collections.defaultdict(list) # {'queryname1':[1, 'ATCGAG', 'TGAGTA', 151], ......}
    bamfile_af = pysam.AlignmentFile(path.realpath(path.expanduser(bamfile)), 'rb')
    query_name_str = ''
    for segment in bamfile_af.fetch():
        try:
            # make the keys of fastq_dict and tags_dict unique
            if segment.is_read1:
                query_name_str = segment.query_name + '/1'
            elif segment.is_read2:
                query_name_str = segment.query_name + '/2'
            else:
                raise NameError

            fastq_dict[query_name_str] = fastq(segment.query_sequence.encode(), segment.query_name, segment.to_dict()['qual'].encode() )
            tags_dict[query_name_str] = [segment.get_tag('Y5'), segment.get_tag('Y3'), segment.get_tag('YL')]

        except Exception as ex:
            common.cprint(segment.query_name)
            print(segment.query_name, ex)
            continue


    # ========================================================correction process begins========================================================
    position_shift_dict = collections.defaultdict(int)  # to record every segment's position shift, it means how many bases we have inserted to a read (a minus value means we deleted some bases)
    last_position_dict = collections.defaultdict(lambda :-1) # to store the last corrected position for a read/segment
    for pileupcolumn in bamfile_af.pileup(max_depth= MAXIMUM_GROUP_MEMBER_AMOUNT, stepper='nofilter', min_base_quality=0 , ignore_overlaps= False):
        if pileupcolumn.nsegments < MINIMUM_GROUP_MEMBER_AMOUNT:
            message = 'Skip position {chrom}:{position}, depth {depth} < {MINIMUM_GROUP_MEMBER_AMOUNT}'.format(chrom= pileupcolumn.reference_name,
                                                                                                             position= pileupcolumn.reference_pos,
                                                                                                             depth= pileupcolumn.nsegments,
                                                                                                             MINIMUM_GROUP_MEMBER_AMOUNT= MINIMUM_GROUP_MEMBER_AMOUNT)
            if VERBOSE:
                print(message)
            #print(pileupcolumn.reference_pos, pileupcolumn.nsegments, pileupcolumn.get_num_aligned() )
            continue

        try:
            message = '@ {chrom}:{pos} Depth:{n}'.format(chrom= pileupcolumn.reference_name,
                                                         pos= pileupcolumn.reference_pos,
                                                         n= pileupcolumn.nsegments)
            #qual_lst = pileupcolumn.get_query_qualities() # a list to store sequencing quality at current position (length equals to coverage)
            #pos_lst = pileupcolumn.get_query_positions() # ditto
            #seq_lst = pileupcolumn.get_query_sequences()   # we have to change the case to upper list
            #name_lst = pileupcolumn.get_query_names() # here name_lst is not end with '/1' or '/2', so it's most likely to contain duplicated items.
            qual_lst = []
            pos_lst = []
            seq_lst = []
            name_lst = []
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    if pileupread.alignment.is_read1:
                        name_lst.append( pileupread.alignment.query_name+'/1' )
                    elif pileupread.alignment.is_read2:
                        name_lst.append( pileupread.alignment.query_name+'/2' )
                    else:
                        raise NameError

                    pos_lst.append( pileupread.query_position )
                    seq_lst.append( pileupread.alignment.query_sequence[pileupread.query_position] )
                    qual_lst.append( pileupread.alignment.query_qualities[pileupread.query_position] )

            #print(len(pileupcolumn.pileups), qual_lst[:20])
            #print('name_lst', len(name_lst), len(set(name_lst)))
            if len(name_lst) != len(set(name_lst)):
                cprint( '========>' + str(len(name_lst)) + ' : ' + str(len(set(name_lst))) )
                print()

            if seq_lst == []: # no coverage at this location
                continue

        except Exception as ex:
            print('error on pileupcolumn: ' + str(ex))
            continue


        #=============================================================================================================================================================================================================
        # deal with insertion here. First determine the consensus insert bases ( a consensus insertion indicates a true insertion happens. If no true insertion happened, consensus insert bases will be an empty string )
        # and then insert the bases between seq[last_pos] and seq[current_pos]
        insert_base_lst = [] # to store insert base for each read/AlignmentSegment, like ['', '', 'A', 'GG', ''], length the equals to coverage at current position
        consensus_insert_b = b'' # the consensus insertion base (DO NOT confuse with consensus base) at current postion. An empty string indicates no true insertion happened.
        for i in range( len(name_lst) ):
            pos_shift_int = position_shift_dict[name_lst[i]]
            insert_base_lst.append( fastq_dict[name_lst[i]].seq[ last_position_dict[name_lst[i]]+pos_shift_int+1 : pos_lst[i]+pos_shift_int ] )  # every item is a bytecode type

        consensus_insert_lst = collections.Counter(insert_base_lst).most_common(2) # [(b'T', 11), (b'A', 10)]
        if len(collections.Counter(insert_base_lst)) == 1:
            consensus_insert_b = b'X'  # we don't make the correction at this position

        elif CORRECTION_VOTING_RATIO is None: # simple majority principle
            consensus_insert_b = consensus_insert_lst[0][0]

        elif isinstance(CORRECTION_VOTING_RATIO, float):

            if consensus_insert_lst[0][1] / len(insert_base_lst) >= CORRECTION_VOTING_RATIO: # if the most frequent one > ration
                consensus_insert_b = consensus_insert_lst[0][0]
            else:
                consensus_insert_b = b'X'  # we don't make the correction at this position

        # make the correction for true/false insertion. Note: the length of consensus_insert_b could be more than 1
        if consensus_insert_b != b'X':
            for i in range( len(name_lst) ):

                pos_shift_int = position_shift_dict[name_lst[i]]
                oriseq_b = fastq_dict[name_lst[i]].seq  # a byte string 'ATCGATCGATCGATCAGCAT....'
                oriqual_ar = fastq_dict[name_lst[i]].qual
                if oriseq_b[last_position_dict[name_lst[i]]+1+pos_shift_int: pos_lst[i]+pos_shift_int] == consensus_insert_b:  # this segment alreay has the correct sequence at this position, skip to next segment
                    continue
                realseq_b = oriseq_b[:last_position_dict[name_lst[i]]+pos_shift_int+1] + consensus_insert_b + oriseq_b[pos_lst[i]+pos_shift_int:] # realseq_b is a true sequence for a read/AlignmentSegment
                realqual_ar = np.concatenate([ oriqual_ar[:last_position_dict[name_lst[i]]+pos_shift_int+1], [qual_lst[i]]*len(consensus_insert_b), oriqual_ar[pos_lst[i]+pos_shift_int:] ]).astype(np.uint8) # generate the true quality array as well
                realqualstr_b = __qualarray2qualstring(realqual_ar)
                fastq_dict[name_lst[i]] = fastq(realseq_b, name_lst[i], realqualstr_b)


                if VERBOSE:
                    message = 'Correct segment {} @ {} to {}'.format(name_lst[i],
                                                                     last_position_dict[name_lst[i]],
                                                                     consensus_insert_b )
                    print(message)
                position_shift_dict[name_lst[i]] += len(consensus_insert_b) - ( pos_lst[i] - last_position_dict[name_lst[i]] - 1) # ( pos[i] - last_position_dict[name_lst[i]] - 1) is origalnal step variation between two pileupcolumns
        #=============================================================================================================================================================================================================

        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # calculate the consensus base (DO NOT confuse with consensus insertion base) at current position
        consensus_lst = collections.Counter( seq_lst ).most_common(2)  # we use a simple majority-consensus method
        if len(consensus_lst) == 1:
            consensus_str = 'X'  # we don't make the correction at this position

        elif CORRECTION_VOTING_RATIO is None: # simple majority principle
            if consensus_lst[0][1] > consensus_lst[1][1]:
                consensus_str = consensus_lst[0][0]
            else: # Two most commons have equal amounts, we calculate and compare the mean sequencing quality/score:
                score_lst = []
                indices = [i for i, x in enumerate(seq_lst) if x == consensus_lst[0][0]]  # which position have top common No. 1 base
                score_lst.append( np.array(qual_lst)[indices].mean() )
                indices = [i for i, x in enumerate(seq_lst) if x == consensus_lst[1][0]]  # which position have top common No. 2 base
                score_lst.append( np.array(qual_lst)[indices].mean() )
                consensus_str = consensus_lst[score_lst[0] < score_lst[1]][0]

        elif isinstance(CORRECTION_VOTING_RATIO, float):
            if consensus_lst[0][1] / len(insert_base_lst) >= CORRECTION_VOTING_RATIO: # if the most frequent one > ration
                if consensus_lst[0][1] > consensus_lst[1][1]:
                    consensus_str = consensus_lst[0][0]
                else: # Two most commons have equal amounts, we calculate and compare the mean sequencing quality/score:
                    score_lst = []
                    indices = [i for i, x in enumerate(seq_lst) if x == consensus_lst[0][0]]  # which position have top common No. 1 base
                    score_lst.append( np.array(qual_lst)[indices].mean() )
                    indices = [i for i, x in enumerate(seq_lst) if x == consensus_lst[1][0]]  # which position have top common No. 2 base
                    score_lst.append( np.array(qual_lst)[indices].mean() )
                    consensus_str = consensus_lst[score_lst[0] < score_lst[1]][0]

            else:
                consensus_str = 'X'

        else:
            pass

        # deal with false insertion, deletion and substitution at mapping position
        for i in range(len(name_lst)):
            if consensus_str != 'X':
            #for i in range(len(name_lst)):
                pos_shift_int = position_shift_dict[name_lst[i]]
                oriseq_b = fastq_dict[name_lst[i]].seq
                oriqual_ar = fastq_dict[name_lst[i]].qual

                if seq_lst[i] == '' and consensus_str != '':                                                 # a false deletion, have to insert the consensus_str and qual-string
                    realseq_b = oriseq_b[:pos_lst[i]+pos_shift_int] + consensus_str.encode() + oriseq_b[pos_lst[i]+pos_shift_int:]
                    #qual_ar = fastq_dict[name_lst[i]].qual
                    realqual_int = int(round(np.mean(qual_lst)))
                    realqual_ar = np.insert( oriqual_ar, pos_lst[i]+pos_shift_int, realqual_int )
                    realqualstr_b = __qualarray2qualstring(realqual_ar)
                    fastq_dict[name_lst[i]] = fastq(realseq_b, name_lst[i], realqualstr_b)
                    if VERBOSE:
                        print('A false deletion in {name} @ {pos}, insert {ins} with quality {qual}'.format(name= name_lst[i],
                                                                                                            pos= pos_lst[i],
                                                                                                            ins= consensus_str,
                                                                                                            qual= realqual_int ))
                    position_shift_dict[name_lst[i]] += 1

                elif seq_lst[i] != '' and consensus_str == '':                                              # a false insertion, have to remove the current seq
                    realseq_b = oriseq_b[:pos_lst[i]+pos_shift_int] + oriseq_b[pos_lst[i]+pos_shift_int+1:]
                    realqual_ar = np.delete(oriqual_ar, pos_lst[i]+pos_shift_int)
                    realqualstr_b = __qualarray2qualstring(realqual_ar)
                    fastq_dict[name_lst[i]] = fastq(realseq_b, name_lst[i], realqualstr_b)
                    if VERBOSE:
                        print('A false insertion in {name} @ {pos}. Removed {rem}'.format( name= name_lst[i],
                                                                                           pos= pos_lst[i],
                                                                                           rem = seq_lst[i] ))
                    position_shift_dict[name_lst[i]] -= 1

                elif consensus_str != '' and seq_lst[i] != '' and seq_lst[i] != consensus_str:              # a substitution/sequencing error
                    realseq_b =  oriseq_b[:pos_lst[i]+pos_shift_int] + consensus_str.encode() + oriseq_b[pos_lst[i]+pos_shift_int+1:]
                    fastq_dict[ name_lst[i] ].seq = realseq_b
                    if VERBOSE:
                        print('A false substitution in {name} @ {pos}. Correct {ori} into {new}'.format(name= name_lst[i],
                                                                                                        pos= pos_lst[i],
                                                                                                        ori= seq_lst[i],
                                                                                                        new= consensus_str))
                    # no need to change the qualstr and position deviation

                else: # reserved got future use
                    pass

            last_position_dict[name_lst[i]] = pos_lst[i]
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    bamfile_af.close()
    with open(output_file, 'a') as results_f:
        read_name_lst = list(fastq_dict.keys())
        read_name_lst.sort()
        if len(read_name_lst) % 2 != 0:
            message = 'read list contains an uneven members.\n' + str(read_name_lst)
            raise IndexError(message)

        for read_name_str in read_name_lst:  # each value is a HTSeq.SequenceWithQualities object
            try:
                if read_name_str[-2:] == '/1':
                    read = fastq_dict[read_name_str]
                elif read_name_str[-2:] == '/2':
                    read = fastq_dict[read_name_str].get_reverse_complement()
                else:
                    message = 'read name in wrong format: {read_name_str}'.format(read_name_str= read_name_str)
                    raise NameError(message)

                read.name = '{read_name_str} Y5:Z:{Y5}\tY3:Z:{Y3}\tYL:i:{YL}\tYA:i:{YA}'.format(read_name_str= read_name_str,
                                                                                                Y5= tags_dict[read_name_str][0],
                                                                                                Y3= tags_dict[read_name_str][1],
                                                                                                YL= tags_dict[read_name_str][2],
                                                                                                YA= 0)
                read.write_to_fastq_file( results_f )

                if VERBOSE:
                    message = 'Write {name}'.format(name = read.name)
                    print(message)

            except Exception as ex:
                cprint(str(ex))
                continue

    return

def deinterleave(fastq: str) -> None:
    '''
    deinterleave a interleaved fastq file
    '''
    forward_f = open(common.add_fix(CORRECTRED_FASTQ, '_1P'), 'w')
    reverse_f = open(common.add_fix(CORRECTRED_FASTQ, '_2P'), 'w')

    for read in HTSeq.FastqReader(fastq):
        if read.name.split()[0][-2:] == '/1':
            read.write_to_fastq_file(forward_f)
        elif read.name.split()[0][-2:] == '/2':
            read.write_to_fastq_file(reverse_f)
        else:
            raise NameError(read.name)

    forward_f.close()
    reverse_f.close()

    return
def main(argvList = sys.argv, argv_int = len(sys.argv)):
    '''
    '''
    t = time.time()
    if readable(SAMPLE_NAME + '_1P'):
        infile_1_str = SAMPLE_NAME + '_1P'
    else:
        message = '{SAMPLE_NAME}_1P not found.'.format(SAMPLE_NAME= SAMPLE_NAME)
        raise FileNotFoundError(message)

    if readable(SAMPLE_NAME + '_2P'):
        infile_2_str = SAMPLE_NAME + '_2P'
    else:
        message = '{SAMPLE_NAME}_2P not found.'.format(SAMPLE_NAME= SAMPLE_NAME)
        raise FileNotFoundError(message)


    #=============================================================
    if not path.isdir('temp'):
        os.mkdir('temp')

    if readable(CORRECTRED_FASTQ):
        os.remove(CORRECTRED_FASTQ)

    if readable('temp/{}.idx'.format(SAMPLE_NAME)):
        os.remove('temp/{}.idx'.format(SAMPLE_NAME))

    #=============preprocess=============

    if __name__ == "__main__":
        p1 = Process(target= FastqPreprocess, args=(infile_1_str, True,) )
        p2 = Process(target= FastqPreprocess, args=(infile_2_str, False, ) )
        print('trimming...')
        p1.start()
        p2.start()
        p1.join()
        p2.join()
        print()

    print('indexing...')

    original_read_dict = SeqIO.index_db('temp/{}.idx'.format(SAMPLE_NAME), filenames=['temp/{}_1P_trimmed'.format(SAMPLE_NAME), 'temp/{}_2P_trimmed'.format(SAMPLE_NAME) ], format='fastq')


    #==============align=============
    forward_file_name = 'temp/' + common.add_fix(path.split(infile_1_str)[1], '_trimmed')
    reverse_file_name = 'temp/' + common.add_fix(path.split(infile_2_str)[1], '_trimmed')
    print('aligning...')

    # do not add '-a' option, we don't want unpaired reads outputed in the results.
    # -U INT, Penalty for an unpaired read pair.
    # -T INT, Don’t output alignment with score lower than INT.
    command_str = 'bwa mem -t {CPU_INT} -C -R "@RG\\tID:{SAMPLE_NAME}\\tSM:{SAMPLE_NAME}\\tPL:illumina" {GENOME_REF} {forward_file_name} {reverse_file_name} > temp/{SAMPLE_NAME}.sam'.format(CPU_INT= CPU_INT,
                                                                                                                                                                                              SAMPLE_NAME= SAMPLE_NAME,
                                                                                                                                                                                              GENOME_REF= GENOME_REF,
                                                                                                                                                                                              forward_file_name= forward_file_name,
                                                                                                                                                                                              reverse_file_name= reverse_file_name)
    if VERBOSE:
        print(command_str)
        print()
    r = common.run(command_str, VERBOSE)
    if r[0] != 0:
        message = 'An error occurs when the following command runs:\n{}'.format(command_str)
        common.cprint(message)
        common.cprint(r[1])
        sys.exit(1)

    # filter and sort by queryname
    # 2316 = 4 + 8 + 256 + 2048 ( UNMAP + MUNMAP + SECONDARY + SUPPLEMENTARY )
    filtered_sam = 'temp/{SAMPLE_NAME}.filtered.sam'.format(SAMPLE_NAME= SAMPLE_NAME)
    with open(filtered_sam, 'w') as f:
        f.write( pysam.view('-h','-f', '1', '-F', '2316', 'temp/{SAMPLE_NAME}.sam'.format(SAMPLE_NAME= SAMPLE_NAME)) )
    queryname_sorted_sam = 'temp/{SAMPLE_NAME}.name.sort.sam'.format(SAMPLE_NAME= SAMPLE_NAME)
    pysam.sort('-n', '-o', queryname_sorted_sam, filtered_sam)


    #=============group, correct and write=============
    queryname_sorted_sam = 'temp/{SAMPLE_NAME}.name.sort.sam'.format(SAMPLE_NAME= SAMPLE_NAME)
    grouped_segments_dict = GroupSegments(queryname_sorted_sam)
    grouped_segments_subdict_list = common.divide_dict(grouped_segments_dict, CPU_INT)

    def __multiprocess_wrapper(grouped_segments_dict:dict, head_dict:dict, part_index:int, counter:multiprocessing.Value, total:int, Lock: multiprocessing.Lock) -> None:
        temp = time.time()
        for group_tag, segment_lst in grouped_segments_dict.items():
            tempfile_str = 'temp/{SAMPLE_NAME}.temp.{part_index}.bam'.format(SAMPLE_NAME= SAMPLE_NAME, part_index= part_index)
            try:
                bamfile_af = pysam.AlignmentFile(tempfile_str, 'wb', header= header_dict)
                for segment in segment_lst:
                    bamfile_af.write(segment)
                bamfile_af.close()
            except:
                message = 'Fail to write segment list {segment_lst}'.format(segment_lst= segment_lst)
                cprint(message)
                continue
            tempsortfile_str = 'temp/{SAMPLE_NAME}.temp.{part_index}.sort.bam'.format(SAMPLE_NAME= SAMPLE_NAME, part_index= part_index)
            pysam.sort('-o', tempsortfile_str, tempfile_str)
            pysam.index(tempsortfile_str)
            tempoutputfile_str = 'temp/{SAMPLE_NAME}.temp.{part_index}.fastq'.format(SAMPLE_NAME= SAMPLE_NAME, part_index= part_index)
            CorrectAndWrite(tempsortfile_str, tempoutputfile_str)
            Lock.acquire()
            counter.value += 1
            message = '[CorrectAndWrite][Elapsed time: {time}s] Corrected groups: {counter}/{total} {group_tag}'.format(time= round(time.time()-temp),
                                                                                                                        counter= counter.value,
                                                                                                                        total= total,
                                                                                                                        group_tag= group_tag)
            print(message, end='\r')
            Lock.release()
        return

    header_dict = pysam.AlignmentFile('temp/{SAMPLE_NAME}.sam'.format(SAMPLE_NAME= SAMPLE_NAME), 'r').header
    counter = Value('i')
    total = len(grouped_segments_dict)
    processed_lst = []
    for part_index in range(CPU_INT):
        tempoutputfile_str = 'temp/{SAMPLE_NAME}.temp.{part_index}.fastq'.format(SAMPLE_NAME= SAMPLE_NAME, part_index= part_index)
        if readable(tempoutputfile_str):
            os.remove(tempoutputfile_str)
        p = Process(target= __multiprocess_wrapper, args=(grouped_segments_subdict_list[part_index], header_dict, part_index, counter, total, Lock()) )
        processed_lst.append(p)
        p.start()

    for p in processed_lst:
        p.join()

    command_str = 'cat '
    for part_index in range(CPU_INT):
        tempoutputfile_str = 'temp/{SAMPLE_NAME}.temp.{part_index}.fastq'.format(SAMPLE_NAME= SAMPLE_NAME, part_index= part_index)
        command_str += tempoutputfile_str + ' '
    command_str += '>> ' + CORRECTRED_FASTQ
    #print(command_str)
    r = common.run(command_str, VERBOSE)
    if r[0] != 0:
        cprint(command_str)
        cprint(r[1] + '\n' + r[2])
        sys.exit(1)

    print('\ndeinterleaving...')
    deinterleave(CORRECTRED_FASTQ)
    os.remove(CORRECTRED_FASTQ)

    print()
    print('done', time.time()-t)
    return

main()
#queryname_sorted_sam = 'temp/{SAMPLE_NAME}.name.sort.sam'.format(SAMPLE_NAME= SAMPLE_NAME)
#grouped_segments_dict = GroupSegments(queryname_sorted_sam)
