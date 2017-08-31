#!/usr/bin/env pyton3

import os
import pandas


#############
# FUNCTIONS #
#############

def parse_key_and_write_config_files(key_file, outdir):
    '''Group key_file rows by 'Flowcell' and 'Lane', and write tsv of 'Barcode'
    and 'Sample' separately for each lane and flowcell to 'outdir' '''
    # generate dicts
    key_data = pandas.read_csv(key_file, delimiter='\t')
    grouped_key_data = key_data.groupby(['Flowcell', 'Lane'])

    # write the key files by group
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    for name, group in grouped_key_data:
        prefix = '_'.join([str(x) for x in name])
        config_file = os.path.join(outdir, '%s.config' % prefix)
        subset = group[['Barcode', 'Sample']]
        if len(subset) > 0:
            subset.to_csv(config_file,
                          sep='\t',
                          header=False,
                          index=False)


###########
# GLOBALS #
###########

key_file = 'data/SQ0003.txt'
outdir = 'output'
stacks_config_dir = os.path.join(outdir, 'stacks_config')

# read key file
key_data = pandas.read_csv(key_file, delimiter='\t')
grouped_key_data = key_data.groupby(['Flowcell', 'Lane'])

# get a list of flowcell_lane
all_fc_lanes = []
for name, group in grouped_key_data:
    all_fc_lanes.append('_'.join([str(x) for x in name]))

#########
# RULES #
#########

# extract per-flowcell/lane sample:barcode information
rule extract_barcode_config:
    input:
        key_file
    output:
        expand('{stacks_config_dir}/{fc_lane}.config',
               stacks_config_dir=stacks_config_dir,
               fc_lane=all_fc_lanes)
    run:
        parse_key_and_write_config_files(key_file, stacks_config_dir)
