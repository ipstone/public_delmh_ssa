#!/usr/bin/env python

# Calculate deletion with microhomology from input TSV file
# This script finds perfect microhomology sequences at deletion breakpoints
import os, argparse
import numpy as np
import pandas as pd
from Bio import SeqIO

def find_mh_1d(direction, seq, start_coord, end_coord, coord2pyidx = -1): 
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
    while (
            left_idx >= 0 and right_idx < len(seq) # search within the bounds of reference sequence
        and 
            abs(left_idx - init_left_idx) <= end_coord - start_coord # one member of homology is within the deletion region
        and 
            seq[left_idx] == seq[right_idx] # homology detected
        ):
        mh += seq[left_idx]
        left_idx += step
        right_idx += step
    if left_idx == init_left_idx:
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
    
    return {'found': found, 'mh': mh, 'left_coords': left_coords, 'right_coords': right_coords}
#

def find_mh_2d(seq, start_coord, end_coord, out_type='dict'):
    # idx is python style
    directions = ['up', 'down']
   
    # mh_dict is from find_mh_1d (including both up/down, all combined as rows)
    mh_dict = {direc: find_mh_1d(direc, seq, start_coord, end_coord) for direc in directions}

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

def find_mh_t2t(del_table, ref_dict): # from input table to output table
    select_cols = ["SAMPLE.TUMOR", "SAMPLE.NORMAL", "variantCaller",
    'CHROM', 'POS', 'REF', 'ALT', 'start_position', 'end_position']

    del_mh = del_table.apply(lambda row: find_mh_2d(
                   ref_dict[row['CHROM']],
                   row['start_position'], 
                   row['end_position'], out_type='srs'), axis=1)
    return pd.concat([del_table[select_cols], del_mh], axis=1)
#

def main():
    parser = argparse.ArgumentParser(description='Calculate microhomology at deletion breakpoints.')
    parser.add_argument('--ref-fasta', required=True,
                        help='Path to reference genome .fasta file')
    parser.add_argument('-i', '--input', required=True,
                        help='Input deletion table in TSV format')
    parser.add_argument('-o', '--output', required=True,
                        help='Output file for microhomology results')
    
    args = parser.parse_args()
    ref_dict = SeqIO.to_dict(SeqIO.parse(args.ref_fasta, 'fasta')) 

    # Load, process del, and output
    del_table = pd.read_table(args.input, low_memory=False, comment='#')
    del_table["del_length"] = del_table.REF.apply(len) - del_table.ALT.apply(len)
    del_table = del_table[ del_table.del_length >0 ]
    del_table["start_position"] = del_table.POS + 1
    del_table["end_position"] = del_table.POS + del_table.del_length

    dm_out = find_mh_t2t(del_table, ref_dict)
    dm_out.to_csv(args.output, index=False, sep='\t')
#

if __name__ == '__main__':
    main()




