#!/usr/bin/env pyton3

import csv
import os
import pandas
import re


#############
# FUNCTIONS #
#############

def parse_key_and_write_config_files(key_file, outdir):
    '''Group key_file rows by 'Flowcell' and 'Lane', and write tsv of 'Barcode'
    and 'Sample' separately for each lane and flowcell to 'outdir' '''
    # generate dicts
    key_data = pandas.read_csv(key_file, delimiter='\t')
    grouped_key_data = key_data.groupby(['Flowcell', 'Lane'])

    # make output directory
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # write the key files by group
    for name, group in grouped_key_data:
        prefix = '_'.join([str(x) for x in name])
        config_file = os.path.join(outdir, '%s.config' % prefix)
        subset = group[['Barcode', 'Sample']]
        if len(subset) > 0:
            subset.to_csv(config_file,
                          sep='\t',
                          header=False,
                          index=False)

    # generate population map
    sample_to_population = {}
    for sample in key_data['Sample']:
        sample_to_population[sample] = re.sub('\d', '', sample).lower()

    # write the population map
    population_map = os.path.join(outdir, 'population_map.txt')
    with open(population_map, 'w') as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows([x, sample_to_population[x]]
                         for x in sample_to_population
                         if x in all_samples)


###########
# GLOBALS #
###########

# files and folders
reference = 'data/genome.fasta'
key_file = 'data/SQ0003.txt'
reads_dir = 'data/raw_reads'
outdir = 'output'
stacks_config_dir = os.path.join(outdir, 'stacks_config')
demultiplex_dir = os.path.join(outdir, 'demux')
population_map = os.path.join(stacks_config_dir, 'population_map.txt')
stacks_dir = os.path.join(outdir, 'stacks')
populations_batch_id = '1'
pop_output = os.path.join(
    stacks_dir,
    'batch_{0}.sumstats.tsv'.format(populations_batch_id))
vcf = os.path.join(
    stacks_dir,
    'batch_{0}.vcf'.format(populations_batch_id))
snprelate_dir = os.path.join(outdir, 'snprelate')
gds = os.path.join(
            snprelate_dir,
            'batch_{}.gds'.format(populations_batch_id))
pca = os.path.join(snprelate_dir,
                   'batch_{}_pca.Rds'.format(populations_batch_id))

stacks_db_dir = os.path.join(outdir, 'stacks_db')
stacks_db_name = "stacks_radtags"

#########
# SETUP #
#########

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
        pca

# extract per-flowcell/lane sample:barcode information
rule extract_barcode_config:
    input:
        key_file
    output:
        expand('{stacks_config_dir}/{fc_lane}.config',
               stacks_config_dir=stacks_config_dir,
               fc_lane=fc_lane_to_sample.keys()),
        population_map
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
            expand('output/demux/{sample}.fq.gz',
                   sample=fc_lane_to_sample[fc_lane]),
        log:
            'output/demux/logs/{0}.log'.format(fc_lane)
        threads:
            1
        shell:
            'bin/stacks/process_radtags '
            '-f {input.read_file} '
            '-i gzfastq -y gzfastq '
            '-b {input.config_file} '
            '-o output/demux '
            '-c -q '
            # '-r --barcode_dist_1 1 '    # rescue barcodes
            '-t 91 '                    # truncate output to 91 b
            '-w 0.1 '                   # window: approx. 9 bases
            '-s 15 '                    # minimum avg PHRED in window
            '--inline_null '
            '--renz_1 apeKI --renz_2 mspI '
            '&> {log}'

# prepare reference genome
rule filter_reference:
    input:
        reference = reference
    output:
        filtered_reference = 'output/genome_index/filtered_reference.fa'
    params:
        min_len = 1000
    threads:
        1
    script:
        'src/filter_reference.py'

rule prepare_reference:
    input:
        'output/genome_index/filtered_reference.fa'
    output:
        'output/genome_index/genome'
    threads:
        1
    log:
        'output/genome_index/index.log'
    shell:
        'bin/gmap/gmap_build '
        '--dir=output/genome_index '
        '--db=genome '
        '{input} '
        '2> {log}'

# map reads per sample
rule map:
    input:
        fq = 'output/demux/{sample}.fq.gz',
        index = 'output/genome_index/genome'
    params:
        genome_dir = 'output/genome_index'
    output:
        'output/map/{sample}.bam'
    threads:
        50
    log:
        'output/map/{sample}.log'
    shell:
        'bin/gmap/gsnap --nthreads={threads} '
        '--npaths=1 '
        '--max-mismatches=5 '
        '--indel-penalty=2 '
        '--min-coverage=0.90 '
        '--format=sam '
        '--db=genome '
        '--dir={params.genome_dir} '
        '{input.fq} '
        '--gunzip '
        '2> {log} '
        '| bin/samtools/samtools sort '
        '-l 9 -m 10G --threads {threads} '
        '--output-fmt BAM '
        '-o {output} '
        '; bin/samtools/samtools index {output}'

# run STACKS
rule stacks:
    input:
        bamfiles = expand('output/map/{sample}.bam',
                          sample=all_samples),
        population_map = population_map
    output:
        pop_output = pop_output,
        vcf = vcf
    params:
        prefix = stacks_dir
    threads:
        50
    log:
        os.path.join(stacks_dir, 'stacks.log')
    run:
        sample_string = ''
        for sample in input.bamfiles:
            sample_string += '-s {} '.format(sample)
        shell('bin/stacks/ref_map.pl -T {threads} '
              '-b {populations_batch_id} '
              '-o {params.prefix} '
              '-O {input.population_map} '
              '-e bin/stacks '
              '{sample_string} '
              '-m 5 '                          # 10 reads per stack
              '-S '                             # disable database
              '-X "populations:-p 6" '         # SNP filtering
              '-X "populations:-r 0.5" '
              '-X "populations:--min_maf 0.1" '
              '-X "populations:--vcf" '        # request populations output
              '-X "populations:--fstats" '
              '-X "populations:--fst_correction p_value" '
              '-X "populations:--kernel_smoothed" '
              '&> {log}')

# re-run populations for FST bootstrapping
rule fst_bootstrap:
    input:
        pop_output = pop_output,
        population_map = population_map
    output:
        'output/populations/batch_1.phistats.tsv'
    params:
        wd = 'output/populations'
    threads:
        50
    log:
        'output/populations/populations.log'
    shell:
        'bin/stacks/populations '
        '-b 1 -P output/stacks '
        '--out_path {params.wd} '
        '-s -t {threads} '
        '-M {input.population_map} '
        '-p 6 -r 0.5 --min_maf 0.1 '
        '--fstats --fst_correction p_value '
        '--kernel_smoothed '
        '--bootstrap '
        '--bootstrap_reps 100 '
        '&> {log}'

rule convert_vcf_to_gds:
    input:
        vcf = vcf
    output:
        gds = gds
    log:
        os.path.join(snprelate_dir,
                     'snprelate.log')
    shell:
        'Rscript -e \''
        'library(SNPRelate) ; '
        'snpgdsVCF2GDS('
        '"{input.vcf}",'
        '"{output.gds}",'
        'method = "biallelic.only",'
        'verbose = TRUE)'
        '\''
        '&> {log}'

rule run_pca:
    input:
        gds = gds
    output:
        pca
    threads:
        50
    log:
        os.path.join(snprelate_dir, 'PCA.log')
    script:
        'src/generate_pca.R'

# after filtering the SNPs properly can easily do a pca with PLINK
# e.g. ~/bin/plink/plink --pca --allow-extra-chr --vcf batch_1.vcf

# create an SQL database for stacks results
rule create_stacks_db:
    input:
        stacks_dir = stacks_dir
    output:
        touch(os.path.join(stacks_db_dir, 'stacks_db.tmp'))
    params:
        db_name = stacks_db_name
    log:
        os.path.join(stacks_db_dir, 'create_db.log')
    shell:
        'mysql -e "CREATE DATABASE {params.db_name}" &> {log} ; '
        'mysql {params.db_name} < bin/stacks/share/stacks/sql/stacks.sql '
        '&>> {log}'

# load stacks results into database (only works on local computer)
rule load_stacks_db:
    input:
        os.path.join(stacks_db_dir, 'stacks_db.tmp')
    output:
        touch(os.path.join(stacks_db_dir, 'stacks_db_load.tmp'))
    params:
        db_name = stacks_db_name,
        stacks_dir = stacks_dir
    log:
        os.path.join(stacks_db_dir, 'load_db.log')
    shell:
        'bin/stacks/bin/load_radtags.pl '
        '-D {params.db_name} '
        '-p {params.stacks_dir} '
        '-b {populations_batch_id} '
        '-M {population_map} '
        '-c '
        '&> {log}'

# index the stacks database for speed
rule index_stacks_db:
    input:
        os.path.join(stacks_db_dir, 'stacks_db_load.tmp')
    output:
        touch(os.path.join(stacks_db_dir, 'stacks_db_index.tmp'))
    params:
        db_name = stacks_db_name
    log:
        os.path.join(stacks_db_dir, 'index_db.log')
    shell:
        'bin/stacks/bin/index_radtags.pl '
        '-s bin/stacks/share/stacks/sql '
        '-D {params.db_name} '
        '-c -t '
        '&> {log}'

# run optimised stacks denovo
rule stacks_denovo:
    input:
        samples = expand('output/demux/{sample}.fq.gz',
                         sample=all_samples),
        sample_dir = 'output/demux',
        population_map = population_map
    output:
        populations_sumstats = ('output/stacks_denovo/'
                                'populations.sumstats.tsv'),
        populations_haplotypes = ('outdir/stacks_denovo/'
                                  'populations.haplotypes.tsv'),
    params:
        wd = 'output/stacks_denovo'
    threads:
        50
    log:
        'output/stacks_denovo/stacks.log'
    shell:
        'bin/stacks/denovo_map.pl '
        '--samples {input.sample_dir} '
        '--popmap {input.population_map} '
        '-T {threads} '
        '-o {params.wd} '
        '-S '
        '-m 3 '
        '-M 3 '
        '-n 3 '
        '-X "populations:-r 0.8"'
        '&> {log} '
