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

##### Run pipeline ###
for sample in samples:
    print "Processing", sample, "....."

    ###########Per-Sample Variant Calling
    #HaplotypeCaller on each sample BAM file 
    #(if a sample's data is spread over more than one BAM, then pass them all in together) to create single-sample gVCFs
    #not recommended for somatic (cancer) variant discovery. For that purpose, use MuTect2 instead
    cmd = ' '.join([gatkCall,  ' -T HaplotypeCaller ', ' -I ' + jp(bamFolder, sample) + '.bam', ' --emitRefConfidence GVCF ', ' -o ' + jp(variantFolder, sample) + '.raw.snps.indels.g.vcf', '>>', logFile, '2>&1'])
    log(cmd, logCommands)
    #os.system(cmd)
    logCommands.close()
    
    '''
    #indel realignment is no longer necessary for variant discovery if you plan to use a variant caller 
    #that performs a haplotype assembly step, such as HaplotypeCaller
    #RealignerTargetCreator
    cmd = ' '.join([gatkCall,  '-T RealignerTargetCreator', ' -I ' + jp(bamFolder, sample) + ".bam", 
                    ' -o ' + jp(bamFolder, sample) + ".RealignerTargetCreator.intervals", '>>', logFile, '2>&1'])
    log(cmd, logCommands)
    #os.system(cmd)
    
    #Realignment Around Indels
    cmd = ' '.join([gatkCall, '-T IndelRealigner', ' -I ' + jp(bamFolder, sample) + ".bam", 
                    ' -o ' + jp(bamFolder, sample) + ".realign.bam", 
                    ' -targetIntervals '  + jp(bamFolder, sample) + ".RealignerTargetCreator.intervals" 
                    '>>', logFile, '2>&1'])
    log(cmd, logCommands)
    #os.system(cmd)

    #Base Quality Score Recalibration
    cmd = ' '.join([gatkCall,  ' -I ' + jp(bamFolder, sample) + ".bam", 
                    ' -o ' + jp(variantFolder, sample) + ".raw.variants.vcf", '>>', logFile, '2>&1'])
    log(cmd, logCommands)
    #os.system(cmd)
'''
'''
for i in samples; do variant=--variant $i;    
    ###########Joint Genotyping
    cmd = ' '.join(['for i in samples; do variant=--variant $i;', gatkCall,  ' -T GenotypeGVCFs ', ' ${variant} ' + jp(variantFolder, sample) + ".raw.snps.indels.g.vcf", 
                    ' -o ' + jp(variantFolder, sample) + "raw.variants.vcf", '>>', logFile, '2>&1'])
    log(cmd, logCommands)
    #os.system(cmd)

    ############Variant Quality Score Recalibration
    #SNPs
    cmd = ' '.join([gatkCall,  ' -T VariantRecalibrator ', ' -I ' + jp(variantFolder, sample) + "raw.variants.vcf",
                    ' -resource:????????????????? '
                    ' -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff ', 
                    ' -mode SNP ', ' -recalFile output.recal ', ' -tranchesFile output.tranches ', ' -rscriptFile output.plots.R ', '>>', logFile, '2>&1'])
    log(cmd, logCommands)
    #os.system(cmd)
    
    #Indels
    cmd = ' '.join([gatkCall,  ' -T VariantRecalibrator ', ' -I ' + jp(variantFolder, sample) + "raw.variants.vcf",
                	' -resource:????????????????? '
                    ' -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff ', 
                    ' -mode INDEL ', ' -recalFile output.recal ', ' -tranchesFile output.tranches ', ' -rscriptFile output.plots.R ', '>>', logFile, '2>&1'])
    log(cmd, logCommands)
    #os.system(cmd)
    
##### Merge all bam files
#cmd = ' '.join(['samtools merge', jp(bamFolder, "*.bam")])
'''