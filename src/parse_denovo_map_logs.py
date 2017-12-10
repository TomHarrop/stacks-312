#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas
import re


def parse_denovo_map_log(dm_log):
    '''
    Parse the dm_log file and return a pandas.DataFrame of sample, unmerged_cov
    and merged_cov 
    '''
    # read the log file
    with open(dm_log, 'rt') as f:
        log_lines = list(enumerate(x.rstrip('\n') for x in f.readlines()))
    # get all the ustacks blocks
    sample_starts = [i for i, line in log_lines
                     if re.match('Sample \d+ of \d+', line)]
    sample_ends = [i for i, line in log_lines
                   if 'ustacks is done' in line]
    # iterate over the blocks and parse the details
    df = pandas.DataFrame()
    for j in range(0, len(sample_starts)):
        sample_lines = list(line for i, line in log_lines
                            if i >= sample_starts[j]
                            and i <= sample_ends[j])
        sample_name = re.sub('Sample \d+ of \d+ \'(?P<sn>.+)\'',
                             '\g<sn>',
                             sample_lines[0])
        unmerged_coverage = [re.sub('.+mean=(?P<mean>\d+\.\d+).+',
                                    '\g<mean>',
                                    line)
                             for line in sample_lines
                             if line.startswith(('Coverage after '
                                                 'assembling stacks:'))][0]
        merged_coverage = [re.sub('.+mean=(?P<mean>\d+\.\d+).+',
                                  '\g<mean>',
                                  line)
                           for line in sample_lines
                           if line.startswith(('Final coverage: '))][0]
        df = df.append({'sample':   sample_name,
                        'unmerged_cov': unmerged_coverage,
                        'merged_cov':   merged_coverage},
                       ignore_index=True)
    # format the data frame
    df = df.reindex(['sample', 'unmerged_cov', 'merged_cov'], axis=1)
    return df


def main():
    # globals
    dm_log = snakemake.input['dm_log']
    stats_output = snakemake.output['coverage_stats']

    # parse the stats
    coverage_stats = parse_denovo_map_log(dm_log)

    # print the results
    coverage_stats.to_csv(
        stats_output,
        index=False)

if __name__ == '__main__':
    main()
