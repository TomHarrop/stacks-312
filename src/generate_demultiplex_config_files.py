#!/usr/bin/env python3

import os
import pandas

# globals
key_file = 'data/SQ0003.txt'
outdir = 'output/stacks_config'


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
