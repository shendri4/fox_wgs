#!/usr/bin/env python
#import argparse
#from glob import glob

from os.path import join as jp
from os.path import abspath
import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-s', "--samples", help="Samples.txt file with sample ID.", required=True)
parser.add_argument('-b', "--bwaindex", help="Path to bwa index file.", required=True)
args = parser.parse_args()
#args = parser.parse_args('-s samples.txt -r /mnt/lfs2/hend6746/fox_cancer/0rawdata_test -b /mnt/lfs2/hend6746/wolves/reference/canfam31/canfam31.fa'.split())

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
        samples.append(l.split('/')[-1].replace('_R1_001.fastq.gz', '').strip())

# Setup folders and paths variables:
bamFolder = abspath('02-Mapped')
variantFolder = abspath('03-Calls')
PBS_scripts = abspath('GATK_PBS_scripts')
#rawdataDir = abspath(args.rawdata)
bwaIndex = abspath(args.bwaindex)
gatkCall = 'java -jar /opt/modules/biology/gatk/3.5/bin/GenomeAnalysisTK.jar -R %s' % bwaIndex

os.system('mkdir -p %s' % bamFolder)
os.system('mkdir -p %s' % variantFolder)
os.system('mkdir -p %s' % PBS_scripts)

##### Run pipeline ###
for sample in samples:
    print "Processing", sample, "....."
    # Set up files:
    logFile = jp(variantFolder, sample + '_GATK.log')
    logCommands = open(jp(PBS_scripts, sample + '_commands.sh'), 'w')

    #Setup for qsub
    log('#!/bin/bash', logCommands)
    log('#PBS -N %s' % sample, logCommands)
    log('#PBS -j oe', logCommands)
    log('#PBS -o %s_job.log' % sample, logCommands)
    log('#PBS -m abe', logCommands)
    log('#PBS -M shendri4@gmail.com', logCommands)
    log('#PBS -q short', logCommands)
    log('#PBS -l mem=100000', logCommands)
    log(". /usr/modules/init/bash", logCommands)
    log("module load python/2.7.10", logCommands)
    log("module load grc", logCommands)

    ###########Per-Sample Variant Calling
    #HaplotypeCaller on each sample BAM file 
    #(if a sample's data is spread over more than one BAM, then pass them all in together) to create single-sample gVCFs
    #not recommended for somatic (cancer) variant discovery. For that purpose, use MuTect2 instead
    cmd = ' '.join([gatkCall,  ' -T HaplotypeCaller ', ' -I ' + jp(bamFolder, sample) + '.bam', ' --emitRefConfidence GVCF ', ' -o ' + jp(variantFolder, sample) + '.raw.snps.indels.g.vcf', '>>', logFile, '2>&1'])
    log(cmd, logCommands)
    #os.system(cmd)
    
for sample in samples:
    variant = '--variant sample' + jp(variantFolder, sample) + '.raw.snps.indels.g.vcf'
    print variant
    ###########Joint Genotyping
    cmd = ' '.join([gatkCall,  ' -T GenotypeGVCFs ', variant, ' -o ' + jp(variantFolder, sample) + 'raw.variants.vcf', '>>', logFile, '2>&1'])
    log(cmd, logCommands)
    
    logCommands.close()
