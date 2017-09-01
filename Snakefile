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
reads_dir = 'data/raw_reads'
outdir = 'output'
stacks_config_dir = os.path.join(outdir, 'stacks_config')
demultiplex_dir = os.path.join(outdir, 'demux')

# read key file
key_data = pandas.read_csv(key_file, delimiter='\t')
grouped_key_data = key_data.groupby(['Flowcell', 'Lane'])

# get a list of fastq files
all_read_files = []
read_dir_files = list((dirpath, filenames)
                      for (dirpath, dirnames, filenames)
                      in os.walk(reads_dir))

for dirpath, filenames in read_dir_files:
    for filename in filenames:
        if 'fastq.gz' in filename:
            all_read_files.append(os.path.join(dirpath, filename))

# get dicts of flowcell_lane:sample and sample:flowcell_lane
fc_lane_to_sample = {}
sample_to_fc_lane = {}
for name, group in grouped_key_data:
    fc_lane = '_'.join([str(x) for x in name])
    sample_list = list(group['Sample'])
    fc_lane_to_sample[fc_lane] = sample_list
    for sample in sample_list:
        sample_to_fc_lane[sample] = fc_lane

# get a list of samples
all_samples = sorted(set(x for x in sample_to_fc_lane.keys()
                         if any(list(sample_to_fc_lane[x] in y
                                     for y in all_read_files))))
all_fc_lanes = [x for x in fc_lane_to_sample
                if any([x in y for y in all_read_files])]

#########
# RULES #
#########

rule all:
    input:
        expand('output/demux/{sample}.fastq.gz',
               sample=all_samples)


# extract per-flowcell/lane sample:barcode information
rule extract_barcode_config:
    input:
        key_file
    output:
        expand('{stacks_config_dir}/{fc_lane}.config',
               stacks_config_dir=stacks_config_dir,
               fc_lane=fc_lane_to_sample.keys())
    run:
        parse_key_and_write_config_files(key_file, stacks_config_dir)

# for loop per fc_lane
for fc_lane in all_fc_lanes:
    rule:
        input:
            read_file = [x for x in all_read_files if fc_lane in x][0],
            config_file = os.path.join(stacks_config_dir,
                                       '%s.config' % fc_lane)
        output:
            expand('output/demux/{sample}.fastq.gz',
                   sample=fc_lane_to_sample[fc_lane]),
            log = 'output/demux/%s/%s.log' % (fc_lane, fc_lane)
        threads:
            1
        shell:
            'bin/stacks/process_radtags '
            '-f {input.read_file} '
            '-i gzfastq -y gzfastq '
            '-b {input.config_file} '
            '-o output/demux '
            '-c -q -r '
            '--inline_null '
            '--renz_1 mspI --renz_2 apeKI'
