#!/usr/bin/env python
#import argparse
#from glob import glob

#-s /mnt/lfs2/hend6746/fox_cancer_test/samples_test.txt
#-r /mnt/lfs2/hend6746/fox_cancer_test/00-RawData
#-b /mnt/lfs2/hend6746/wolves/reference/canfam31/canfam31.fa

#Preparation for reference file
#bwa index -a bwtsw sarHar1.fa (takes awhile)
#samtools faidx sarHar1.fa
#java -jar /mnt/lfs2/hend6746/modules/picard-tools/1.115/CreateSequenceDictionary.jar REFERENCE=sarHar1.fa OUTPUT=sarHar1.dict 

from os.path import join as jp
from os.path import abspath
import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-s', "--samples", help="Samples.txt file with sample ID.", required=True)
parser.add_argument('-r', "--rawdata", help="Path to raw fastq data.", required=True)
parser.add_argument('-b', "--bwaindex", help="Path to bwa index file.", required=True)
args = parser.parse_args()

VERBOSE=False

#Function definitions:
def log(txt, out):
    if VERBOSE:
        print(txt)
    out.write(txt+'\n')
    out.flush()

## Read in samples and put them in a list:
samples = []
for l in open(args.samples):
    if len(l) > 1:
        samples.append(l.split('/')[-1].replace('_L001_R1_001.fastq.gz', '').strip())

# Setup folders and paths variables:
resultsDir = abspath('01-Cleaned')
bamFolder = abspath('02-Mapped')
variantFolder = abspath('03-Calls')
PBS_scripts = abspath('A-cleaning_scripts')
rawdataDir = abspath(args.rawdata)
bwaIndex = abspath(args.bwaindex)
picardCall = 'java -Xms4g -jar /mnt/lfs2/hend6746/modules/picard-tools/1.115/MarkDuplicates.jar '

os.system('mkdir -p %s' % resultsDir)
os.system('mkdir -p %s' % bamFolder)
os.system('mkdir -p %s' % variantFolder)
os.system('mkdir -p %s' % PBS_scripts)

##### Run pipeline ###
for sample in samples:
    print "Processing", sample, "....."
    # Set up files:
    logFile = jp(resultsDir, sample + '_cleaning.log')
    logCommands = open(jp(PBS_scripts, sample + '_cleaning_commands.sh'), 'w')

    #Setup for qsub
    log('#!/bin/bash', logCommands)
    log('#PBS -N %s' % sample, logCommands)
    log('#PBS -j oe', logCommands)
    log('#PBS -o %s_job.log' % sample, logCommands)
    log('#PBS -m abe', logCommands)
    log('#PBS -M shendri4@gmail.com', logCommands)
    log('#PBS -q reg', logCommands)
    log('#PBS -l mem=100gb', logCommands)
    log(". /usr/modules/init/bash", logCommands)
    log("module load python/2.7.10", logCommands)
    log("module load bwa", logCommands)
    log("module load grc", logCommands)

# Second run flash2
# the --max-overlap was set to 600, but that seems really long; default is 65, try 300
# --max-overlap 400 --min-overlap 15 --max-mismatch-density .10 --min-overlap-outie 35 --percent-cutoff 25
    cmd = ' '.join(['flash2 --max-overlap 300 --allow-outies --threads 7', ' -d ', resultsDir, ' -o ', jp(sample + '_flash'),
                    jp(rawdataDir, sample + '.fastq.1.gz'), jp(rawdataDir, sample + '.fastq.2.gz'),
                    '>>', logFile, '2>&1'])
    log(cmd, logCommands)

# Third run sickle
# The --length-threshold was set to 200, but that seems really long; default is 20
    cmd = ' '.join(['sickle pe --length-threshold 20 --qual-threshold 25 --qual-type sanger -f', jp(resultsDir, sample + '_flash.notCombined_1.fastq'),
                    '-r', jp(resultsDir, sample + '_flash.notCombined_2.fastq'),
                    '--output-pe1', jp(resultsDir, sample + '_sickle_PE1.fastq'),
                    '--output-pe2', jp(resultsDir, sample + '_sickle_PE2.fastq'),
                    '--output-single', jp(resultsDir, sample + '_sickle_SE.fastq'), '>>', logFile, '2>&1'])
    log(cmd, logCommands)# 

# Combine SE files:
    cmd = ' '.join(['cat', jp(resultsDir, sample + '_sickle_SE.fastq'), jp(resultsDir, sample + '_flash.extendedFrags.fastq'),
                    '>', jp(resultsDir, sample + "_cleaned_SE.fastq")])
    log(cmd, logCommands)

# Rename PE and SE files to something nicer:
    cmd = ' '.join(['mv', jp(resultsDir, sample + "_sickle_PE1.fastq"), jp(resultsDir, sample + "_cleaned_PE1.fastq")])
    log(cmd, logCommands)

    cmd = ' '.join(['mv', jp(resultsDir, sample + "_sickle_PE2.fastq"), jp(resultsDir, sample + "_cleaned_PE2.fastq")])
    log(cmd, logCommands)

# Compress cleaned files:
    cmd = ' '.join(['gzip', jp(resultsDir, sample + '*.fastq')])
    log(cmd, logCommands)

#     Run BWA to map samples, combine sam files, sort
#     logFile = jp(bamFolder, sample + '_mapping.log')
#     cmd = ' '.join(["bwa mem -t 16 -R '@RG\tID:bwa\tSM:" + sample + "\tPL:ILLUMINA'",
#                     bwaIndex, jp(resultsDir, sample + "_cleaned_PE1.fastq.gz"),
#                     jp(resultsDir, sample + "_cleaned_PE2.fastq.gz"), ">", jp(bamFolder, sample + "_PE.sam"),
#                     "2>", logFile])
#     log(cmd, logCommands)

# Run BWA to map PE samples to reference genome (make sure ref is properly formated (see top of script))
# -t number of threads -R read group header
    logFile = jp(resultsDir, sample + '_mapping.log')
    cmd = ' '.join(["bwa mem -t 16 -R '@RG\tID:bwa\tSM:" + sample + "\tPL:ILLUMINA'",
                    bwaIndex, jp(resultsDir, sample + "_cleaned_PE1.fastq.gz"),
                    jp(resultsDir, sample + "_cleaned_PE2.fastq.gz"), "| samtools view -bS -@ 30 -o", jp(bamFolder, sample) + "PE.bam",
                    "2>", logFile])
    log(cmd, logCommands)

# Run BWA to map SE samples to reference genome (make sure ref is properly formated (see top of script))
    cmd = ' '.join(["bwa mem -t 16 -R '@RG\tID:bwa\tSM:" + sample + "\tPL:ILLUMINA'",
                    bwaIndex, jp(resultsDir, sample + "_cleaned_SE.fastq.gz"), "| samtools view -bS -@ 30 -o", jp(bamFolder, sample + "_SE.bam"),
                    "2>>", logFile])
    log(cmd, logCommands)

# merge
    cmd = ' '.join(['samtools merge -c', jp(bamFolder, sample + ".bam"), jp(bamFolder, sample + "PE.bam"), jp(bamFolder, sample + "SE.bam")])
    log(cmd, logCommands)

# make sure there can be lots of files or it will not be able to handle the samtools sort
    cmd =' '.join(['ulimit -n 2048'])
    log(cmd, logCommands)
    
# sort bam file; -@ number of threads
    cmd = ' '.join(['samtools sort -o', jp(bamFolder, sample) + "sorted.bam", ' -@ 30', jp(bamFolder, sample + ".bam")])
    log(cmd, logCommands)

# Mark PCR duplicates (remove duplicates, if desired)
    cmd = ' '.join([picardCall, ' INPUT=' + jp(bamFolder, sample + "sorted.bam"), ' OUTPUT=' + jp(bamFolder, sample + "_markdup.bam"),
                    ' METRICS_FILE=' + jp(bamFolder, sample + ".metrics"), ' REMOVE_DUPLICATES=true ',
                    ' ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT', '>>', logFile, '2>&1'])
    log(cmd, logCommands)

# Index bam file:
    cmd = ' '.join(['samtools index', jp(bamFolder, sample) + "_markdup.bam"])
    log(cmd, logCommands)     
    logCommands.close()
