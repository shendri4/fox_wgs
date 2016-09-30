#!/usr/bin/env python
#import argparse
#from glob import glob

#Usage: python ../fox_wgs/jointgenotyping_GATK.py -s ./samples.txt -b /mnt/lfs2/hend6746/wolves/reference/canfam31/canfam31.fa

from os.path import join as jp
from os.path import abspath
import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-s', "--samples", help="Samples.txt file with sample ID.", required=True)
parser.add_argument('-b', "--bwaindex", help="Path to bwa index file.", required=True)
args = parser.parse_args()
#args = parser.parse_args('-s samples.txt -r /mnt/lfs2/hend6746/fox_cancer/testdata -b /mnt/lfs2/hend6746/wolves/reference/canfam31/canfam31.fa'.split())

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
PBS_scripts = abspath('jointGATK_PBS_scripts')
#rawdataDir = abspath(args.rawdata)
bwaIndex = abspath(args.bwaindex)
gatkCall = 'java -jar /opt/modules/biology/gatk/3.5/bin/GenomeAnalysisTK.jar -R %s' % bwaIndex

os.system('mkdir -p %s' % bamFolder)
os.system('mkdir -p %s' % variantFolder)
os.system('mkdir -p %s' % PBS_scripts)

logFile = jp(variantFolder, 'joint_GATK.log')
logCommands = open(jp(PBS_scripts, 'joint_commands.sh'), 'w')

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

variants = []
for l in samples:
	if len(l) > 1:
		variants.append(l.join(['--variant' + jp(variantFolder, sample) + '.raw.snps.indels.g.vcf'].strip('/n').split('\t'))
#		line=line.strip('\n')
#		lineList=line.split('\t')

print samples
#+ jp(variantFolder, sample) + '.raw.snps.indels.g.vcf'
###########Joint Genotyping
cmd = ' '.join([gatkCall,  ' -T GenotypeGVCFs ', variant, ' -o ' + jp(variantFolder, sample) + 'raw.variants.vcf', '>>', logFile, '2>&1'])
log(cmd, logCommands)
#os.system(cmd)
logCommands.close()
