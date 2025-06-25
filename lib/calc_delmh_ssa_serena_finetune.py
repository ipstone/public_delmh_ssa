#!/usr/bin/env python
# This script is especially tailored for the Serena's deletion format where the deleted sequence is given as reference_genome_allele (or mutated_from_allele) and mutated to allele is noted as '-' 
# -- we are further fine tunning the definition so that the last nucleotide in
# the matched sequence must be ATGC, not '-'


# Calculate deletion with homeology (imperfect homology) using SSA analysis
# This script finds the longest homeology sequences at deletion breakpoints
import os, argparse
import numpy as np
import pandas as pd
from Bio import SeqIO

# Specify the level of homology allowed
# HOMOLOGY_CUTOFF = 0.8
def find_mh_1d(direction, seq, start_coord, end_coord, coord2pyidx = -1, HOMOLOGY_CUTOFF = 0.8): 
    """arguments: 
        seq - is the chromosome seq object
        start_coord: pos of the first deleted nt (1 indexed)
        end_coord: pos of the last deleted nt (1 indexed)
        as python is 0 indexed, conversion would be:
               python_index = coord - 1; coord = python_index + 1
    """
    if direction == 'up':
        step = -1
        offset = 0
    elif direction == 'down':
        step = 1
        offset = 1
    else:
        raise Exception('`direction` must be either \'up\' or \'down\'.')
    
    mh = ''
    left_idx = init_left_idx = int(start_coord) + coord2pyidx - 1 + offset # numpy.int64 is "Invalid index"
    right_idx = init_right_idx = int(end_coord) + coord2pyidx + offset

    # Previous code finding the exact max match
    # while (
    #         left_idx >= 0 and right_idx < len(seq) # search within the bounds of reference sequence
    #     and 
    #         abs(left_idx - init_left_idx) <= end_coord - start_coord # one member of homology is within the deletion region
    #     and 
    #         seq[left_idx] == seq[right_idx] # homology detected
    #     ):
    #     mh += seq[left_idx]
    #     left_idx += step
    #     right_idx += step

    # Isaac new implementation of homology search for longest match
    for i in range((abs(end_coord-start_coord)+1)):
        if seq[left_idx] == seq[right_idx]: # homology detected
            mh += seq[left_idx]
        else:
            mh +='-'  # Use '-' to indicate no homology

        left_idx +=step
        right_idx +=step

    # Going through mh string, to find the longest sub-string that with homology > HOMOLOGY_CUTOFF
    
    # Find the longest sub-string with homology > HOMOLOGY_CUTOFF
    max_length_homology = 0
    max_length_homology_string = ""
    max_length_homology_start = 0

    for i in range(len(mh)):
        a_string=mh[:(i+1)]
        a_string_homology =  a_string.replace('-','')
        homology_ratio = len(a_string_homology) / len(a_string)
        homology_total_length = len(a_string)

        # Keep the max length homology string
        if ( homology_ratio >= HOMOLOGY_CUTOFF ) and ( homology_total_length > max_length_homology ) and ( a_string[-1] != '-' ):  # last nt must be ATGC
            max_length_homology = homology_total_length
            max_length_homology_string = a_string
            max_length_homology_start = i 

    # Update the left and right index to the start of the max length homology string    
    left_idx = init_left_idx + max_length_homology_start*step
    right_idx = init_right_idx + max_length_homology_start*step

    # Code logic to report the index and max length homology string
    if max_length_homology == 0:
        found = False
        left_coords = right_coords = []
    else:
        found = True
        coords = [idx - coord2pyidx for idx in [init_left_idx, left_idx - step, init_right_idx, right_idx - step]]
        left_coords = coords[:2]
        right_coords = coords[2:]
        if direction == 'up':
            mh = mh[::-1]
            left_coords = left_coords[::-1]
            right_coords = right_coords[::-1]

    # The returned max_length_homology_string is the longest homology string within the homology cutoff
    return {'found': found, 'mh': max_length_homology_string, 'left_coords': left_coords, 'right_coords': right_coords}

def find_mh_2d(seq, start_coord, end_coord, out_type='dict', homology_cutoff=0.8):
    # idx is python style
    directions = ['up', 'down']
   
    # mh_dict is from find_mh_1d (including both up/down, all combined as rows)
    mh_dict = {direc: find_mh_1d(direc, seq, start_coord, end_coord, HOMOLOGY_CUTOFF=homology_cutoff) for direc in directions}

    # mh_srs is calculating up/down, adding results as columns
    mh_list = []
    for direc in directions:
        mh_dict_direc = mh_dict[direc]
        mh_srs_dict_direc = {key: mh_dict_direc[key] for key in ['found', 'mh']}
        for position in ['left', 'right']:
            coords = mh_dict_direc[position + '_coords']
            if len(coords) == 0:
                coords = [np.nan] * 2
            pos_dict = {position + '_start': coords[0], position + '_end': coords[1]}
            mh_srs_dict_direc.update(pos_dict)
        mh_srs_direc = pd.Series(mh_srs_dict_direc)
        mh_srs_direc.index = ['_'.join([direc, idx]) for idx in mh_srs_direc.index]
        mh_list.append(mh_srs_direc)

    mh_srs = pd.concat(mh_list, axis=0) # srs is for both directions
    
    if out_type == 'dict': # result as multple rows
        return mh_dict
    elif out_type == 'srs': # up/down result as multiple columns
        return mh_srs
    else:
        raise Exception('`out_type` must be either \'dict\' or \'srs\'.')
#

def find_mh_t2t(del_table, ref_dict, homology_cutoff): # from input table to output table
    select_cols = ["SAMPLE.TUMOR", "SAMPLE.NORMAL", "variantCaller",
    'CHROM', 'POS', 'REF', 'ALT', 'start_position', 'end_position']

    del_mh = del_table.apply(lambda row: find_mh_2d(
                   ref_dict[row['CHROM']],
                   row['start_position'], 
                   row['end_position'], out_type='srs',
                   homology_cutoff=  homology_cutoff), axis=1)
    return pd.concat([del_table[select_cols], del_mh], axis=1)
#

def main():
    # Configure the genome reference and homology cutoff as arguments
    parser = argparse.ArgumentParser(description='Calculate microhomology.')
    parser.add_argument('--ref-fasta', 
                        help='Path to reference genome .fasta. Must be specified, no default.')

    # Add homology_cutoff as an argument
    parser.add_argument('--homology-cutoff', 
                        help='Minimum homology length to report. Default is 0.8.', 
                        default=0.8, type=float)
    # Adding input deletion table argument 
    parser.add_argument('-i', '--input', required=True,
                        help='Input deletion table in TSV format')
    # Adding output argument
    parser.add_argument('-o', '--output', required=True,
                        help='Output file for homeology analysis results')

    args = parser.parse_args()
    ref_dict = SeqIO.to_dict(SeqIO.parse(args.ref_fasta, 'fasta')) 

    # Load, process del, and output
    del_table = pd.read_table(args.input, low_memory=False, comment='#')
    del_table["del_length"] = del_table.REF.apply(len) - del_table.ALT.apply(len) + 1
    del_table = del_table[ del_table.del_length >0 ]

    # The start_position is different from jrflab output, as the start position is the first position (nucleotide) that's deleted.
    # -- in jrflab modules, the POS, the start_position is the position before the deletion
    del_table["start_position"] = del_table.POS 
    del_table["end_position"] = del_table.POS + del_table.del_length

    dm_out = find_mh_t2t(del_table, ref_dict, homology_cutoff=args.homology_cutoff)
    dm_out.to_csv(args.output, index=False, sep='\t')
#

if __name__ == '__main__':
    main()




