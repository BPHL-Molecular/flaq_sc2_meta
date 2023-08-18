#!/usr/bin/env python


'''
This program takes in Illumina paired-end fastqs using the ARTIC primer schemes for SARS-CoV-2
and reports out quality metrics and variants from metagenomic clinical or environmental SARS-CoV-2 samples
including wastewater samples.
'''

import os
import sys
import subprocess
import argparse
import datetime
import pandas as pd
import re
import os.path

#Parse arguments, get path for fastqs, primer version
parser = argparse.ArgumentParser(usage='flaq_sc2_meta.py <input_dir> [options]')
parser.add_argument('input', help='path to dir with fastqs')
parser.add_argument('--primer_bed', help='path to ARTIC SC2 primer bed file')
parser.add_argument('--lib_frag', default='no_frag', choices=['no_frag', 'frag'], help='specify if input amplicons were fragmented, (default: no_frag)') 
parser.add_argument('--threads', default=8, dest='threads', help='specify number of threads, (default: %(default)s)')
parser.add_argument('--ref_fasta', help='path to reference fasta')
parser.add_argument('--ref_gff', help='path to reference gff')


if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

args = parser.parse_args()

input_dir = os.path.abspath(args.input) + '/'
primers = os.path.abspath(args.primer_bed)
frag = args.lib_frag
threads = str(args.threads)
ref = os.path.abspath(args.ref_fasta)
gff = os.path.abspath(args.ref_gff)
cwd = os.getcwd() + '/'

output_dir = cwd + datetime.date.today().strftime('%Y-%m-%d') + '_flaq_ww_run'
subprocess.run('mkdir -p ' + output_dir, shell=True, check=True) #make output directory date_flaq_legion_run

adapters = '/bbmap/resources/adapters.fa'
phix = '/bbmap/resources/phix174_ill.ref.fa.gz'

#Get sample names
samples = []
fastqs = []

#Look at some code examples to get fastq names R1, _1 or R1_001 (make work for more sample types later)
for f in os.listdir(input_dir):
    if f.endswith('.fastq.gz'):
        fastqs.append(f)
        sn = f.split("_")
        sn = sn[0]
        samples.append(sn)
unique = set(samples)
samples = list(unique)
samples.sort()

#Create output file
report = open(output_dir + '/report.txt', 'w')
header = ['sampleID', 'reference', 'start', 'end', 'num_raw_reads', 'num_clean_reads', 'num_mapped_reads', 'cov_bases_mapped', 'percent_genome_cov_map', 'mean_depth', 'mean_base_qual', 'mean_map_qual', 'freyja_summary', 'freyja_lineage', 'freya_lineage_abund']
#header = ['sampleID', 'reference', 'start', 'end', 'num_raw_reads', 'num_clean_reads', 'num_mapped_reads', 'cov_bases_mapped', 'percent_genome_cov_map', 'mean_depth', 'mean_base_qual', 'mean_map_qual']
report.write('\t'.join(map(str,header)) + '\n')

#Run pipeline for each sample
for s in samples:
    sample_dir = output_dir + '/' + s + '/'
    subprocess.run('mkdir -p ' + sample_dir, shell=True, check=True) #mkdir for each sample name
    subprocess.run('cp ' + input_dir + s + '*.fastq.gz ' + sample_dir, shell=True, check=True) #cp fastqs to new dir

    out_log = open(sample_dir + s + '.out', 'w')
    err_log = open(sample_dir + s + '.err', 'w')

    #Get number of raw reads
    proc_1 = subprocess.run('zcat ' + sample_dir + s + '_1.fastq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
    wc_out_1 = proc_1.stdout.rstrip()
    reads_1 = int(wc_out_1) / 4
    proc_2 = subprocess.run('zcat ' + sample_dir + s + '_2.fastq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
    wc_out_2 = proc_2.stdout.rstrip()
    reads_2 = int(wc_out_2) / 4
    raw_reads = reads_1 + reads_2
    raw_reads = int(raw_reads)

    #Run fastqc on original reads 
    subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/fastqc_0.11.9.sif fastqc ' + sample_dir + '*.fastq.gz --threads ' + threads + ' --outdir ' + sample_dir, shell=True, stdout=out_log, stderr=err_log, check=True)
    #Rename fastqc output files
    subprocess.run('mv ' + sample_dir + s + '_1_fastqc.html ' + sample_dir + s + '_1_original_fastqc.html', shell=True, check=True)
    subprocess.run('mv ' + sample_dir + s + '_1_fastqc.zip ' + sample_dir + s + '_1_original_fastqc.zip', shell=True, check=True)
    subprocess.run('mv ' + sample_dir + s + '_2_fastqc.html ' + sample_dir + s + '_2_original_fastqc.html', shell=True, check=True)
    subprocess.run('mv ' + sample_dir + s + '_2_fastqc.zip ' + sample_dir + s + '_2_original_fastqc.zip', shell=True, check=True)

    #Run trimmomatic
    subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/trimmomatic_0.39.sif trimmomatic PE -threads ' + threads + ' -phred33 -trimlog ' + sample_dir + s + '.log ' + sample_dir + s + '_1.fastq.gz ' + sample_dir + s + '_2.fastq.gz ' + sample_dir + s + '_1.trim.fq.gz ' + sample_dir + s + '_unpaired_1.trim.fq.gz ' + sample_dir + s + '_2.trim.fq.gz ' + sample_dir + s + '_unpaired_2.trim.fq.gz SLIDINGWINDOW:4:30 MINLEN:75 TRAILING:20 > ' + sample_dir + s + '_trimstats.txt', shell=True, stdout=out_log, stderr=err_log, check=True)
    #rm unpaired reads
    subprocess.run('rm ' + sample_dir + s + '*_unpaired*.trim.fq.gz', shell=True, check=True)
    #rm fastq files copied from previous dir
    subprocess.run('rm ' + sample_dir + s + '*.fastq.gz', shell=True, check=True)

    #Run bbduk to remove Illumina adapter sequences and any PhiX contamination
    subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/bbtools_38.94.sif bbduk.sh in1=' + sample_dir + s + '_1.trim.fq.gz in2=' + sample_dir + s + '_2.trim.fq.gz out1=' + sample_dir + s + '_1.rmadpt.fq.gz out2=' +sample_dir + s + '_2.rmadpt.fq.gz ref=' + adapters + ' stats=' + sample_dir + s + '.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo', shell=True, stdout=out_log, stderr=err_log, check=True)
    subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/bbtools_38.94.sif bbduk.sh in1=' + sample_dir + s + '_1.rmadpt.fq.gz in2=' + sample_dir + s + '_2.rmadpt.fq.gz out1=' + sample_dir + s + '_1.fq.gz out2=' + sample_dir + s + '_2.fq.gz outm=' + sample_dir + s + '_matchedphix.fq ref=' + phix + ' k=31 hdist=1 stats=' + sample_dir + s + '_phixstats.txt', shell=True, stdout=out_log, stderr=err_log, check=True)
    subprocess.run('rm ' + sample_dir + '*.trim.fq.gz', shell=True, check=True)
    subprocess.run('rm ' + sample_dir + '*.rmadpt.fq.gz', shell=True, check=True)

    #Run fastqc on clean forward and reverse reads
    subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/fastqc_0.11.9.sif fastqc ' + sample_dir + '*.fq.gz --threads ' + threads, shell=True, stdout=out_log, stderr=err_log, check=True)
    #Rename fastqc output files
    subprocess.run('mv ' + sample_dir + s + '_1_fastqc.html ' + sample_dir + s + '_1_clean_fastqc.html', shell=True, check=True)
    subprocess.run('mv ' + sample_dir + s + '_1_fastqc.zip ' + sample_dir + s + '_1_clean_fastqc.zip', shell=True, check=True)
    subprocess.run('mv ' + sample_dir + s + '_2_fastqc.html ' + sample_dir + s + '_2_clean_fastqc.html', shell=True, check=True)
    subprocess.run('mv ' + sample_dir + s + '_2_fastqc.zip ' + sample_dir + s + '_2_clean_fastqc.zip', shell=True, check=True)

    #Run multiqc
    subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/multiqc_1.8.sif multiqc ' + sample_dir + '*_fastqc.zip -o ' + sample_dir, shell=True, stdout=out_log, stderr=err_log, check=True)

    #Get number of clean reads
    proc_c1 = subprocess.run('zcat ' + sample_dir + s + '_1.fq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
    wc_out_c1 = proc_c1.stdout.rstrip()
    reads_c1 = int(wc_out_c1) / 4
    proc_c2 = subprocess.run('zcat ' + sample_dir + s + '_2.fq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
    wc_out_c2 = proc_c2.stdout.rstrip()
    reads_c2 = int(wc_out_c2) / 4
    clean_reads = reads_c1 + reads_c2
    clean_reads = int(clean_reads)

    #Map reads to reference
    align_dir = sample_dir + 'alignment/'
    subprocess.run('mkdir ' + align_dir, shell=True, check=True)

    #If frag == 'no_frag', do not remove PCR duplicates
    if frag == 'no_frag':
        subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/bwa_0.7.17.sif bwa mem -t ' + threads + ' ' + ref + ' ' + sample_dir + s + '_1.fq.gz ' + sample_dir + s + '_2.fq.gz | singularity exe\
c -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools view - -F 4 -u -h | singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools sort > ' + align_dir + s + '.sorted.bam', shell=True, check=True)
    #If frag == 'frag', remove PCR duplicates
    elif frag == 'frag':
        subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/bwa_0.7.17.sif bwa mem -t ' + threads + ' ' + ref + ' ' + sample_dir + s + '_1.fq.gz ' + sample_dir + s + '_2.fq.gz | singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools view - -F 4 -u -h | singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools sort -n > ' + align_dir + s + '.namesorted.bam', shell=True, check=True) #output name sorted bam
        subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools fixmate -m ' + align_dir + s + '.namesorted.bam ' + align_dir + s + '.fixmate.bam', shell=True, check=True) #fixmate
        #Create positional sorted bam from fixmate.bam
        subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools sort -o ' + align_dir + s + '.positionsort.bam ' + align_dir + s + '.fixmate.bam', shell=True, check=True)
        #Mark duplicate reads
        subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools markdup ' + align_dir + s + '.positionsort.bam ' + align_dir + s + '.markdup.bam', shell=True, check=True)
        #Remove duplicate reads
        subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools markdup -r ' + align_dir + s + '.positionsort.bam ' + align_dir + s + '.dedup.bam', shell=True, check=True)
        #Sort dedup.bam and rename to .sorted.bam
        subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools sort -o ' + align_dir + s + '.sorted.bam ' + align_dir + s + '.dedup.bam', shell=True, check=True)


    #Index final sorted bam from either no_frag or frag paths
    subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools index ' + align_dir + s + '.sorted.bam', shell=True, check=True)

    #Trim primers with iVar
    subprocess.run('cd ' + align_dir + ' && singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/ivar_1.3.1.sif ivar trim -i ' + s + '.sorted.bam  -b ' + primers + ' -p ' + s + '.primertrim -e', shell=True, stdout=out_log, stderr=err_log, check=True)
    subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools sort ' + align_dir + s + '.primertrim.bam -o ' + align_dir + s + '.primertrim.sorted.bam', shell=True, stdout=out_log, stderr=err_log, check=True)
    subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools coverage ' + align_dir  + s + '.primertrim.sorted.bam -o ' + align_dir + s + '.coverage.txt', shell=True, stdout=out_log, stderr=err_log, check=True)

    #Get map stats
    with open(align_dir + s + '.coverage.txt', 'r') as cov_report:
        header = cov_report.readline()
        header = header.rstrip()
        stats = cov_report.readline()
        stats = stats.rstrip()
        stats = stats.split()
        ref_name = stats[0]
        start = stats[1]
        end = stats[2]
        reads_mapped = stats[3]
        cov_bases = stats[4]
        cov = stats[5]
        depth = stats[6]
        baseq = stats[7]
        mapq = stats[8]

    #Call variants
    var_path = sample_dir + 'variants/'
    subprocess.run('mkdir ' + var_path, shell=True, check=True)
    subprocess.run('cd ' + var_path + ' && singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools mpileup -A -d 1000000 --reference ' + ref + ' -B -Q 0 ' + align_dir + s + '.primertrim.sorted.bam | singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/ivar_1.3.1.sif ivar variants -r ' + ref + ' -m 10 -p ' + s + '.variants -q 20 -t 0.03 -g ' + gff, shell=True, stdout=out_log, stderr=err_log, check=True)

    #Run Freyja
    subprocess.run('freyja update', shell=True, check=True, stdout=out_log, stderr=err_log)
    subprocess.run('mkdir ' + align_dir + 'freyja/', shell=True, check=True)
    subprocess.run('freyja variants ' + align_dir + s + '.primertrim.sorted.bam --variants ' + align_dir + 'freyja/' + s + '.variants --depths ' + align_dir + 'freyja/' + s + '.depths --ref ' + ref, shell=True, stdout=out_log, stderr=err_log, check=True)
    subprocess.run('freyja demix --eps 0.01 ' + align_dir + 'freyja/' + s + '.variants.tsv ' + align_dir + 'freyja/' + s + '.depths --output ' + align_dir + 'freyja/' + s + '.freyja.out', shell=True, stdout=out_log, stderr=err_log, check=True)

# ADD BACK IN ONCE GITHUB ISSUE RESOLVED FOR OUTPUT FORMAT
#    #Parse freyja
    with open(align_dir + 'freyja/' + s + '.freyja.out', 'r') as freyja_out:
        header = freyja_out.readline()
        summ = freyja_out.readline()
        summ = summ.rstrip()
        summ = summ.split("\t")[1]
        lineage = freyja_out.readline()
        lineage = lineage.rstrip()
        lineage = lineage.split("\t")[1]
        abund = freyja_out.readline()
        abund = abund.rstrip()
        abund = abund.split("\t")[1]

    #Write to output file
    results = [s, ref_name, start, end, raw_reads, clean_reads, reads_mapped, cov_bases, cov, depth, baseq, mapq, summ, lineage, abund] 
#    results = [s, ref_name, start, end, raw_reads, clean_reads, reads_mapped, cov_bases, cov, depth, baseq, mapq]
    report.write('\t'.join(map(str,results)) + '\n')
    out_log.close()
    err_log.close()

 
report.close()
