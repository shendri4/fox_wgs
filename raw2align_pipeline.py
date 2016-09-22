#!/usr/bin/env python
#import argparse
#from glob import glob
'''
Make sure to load:
module load python/2.7.10
module load bwa
module load grc
module load samtools
'''

from os.path import join as jp
import os

VERBOSE=False

#Function definitions:
def log(txt, out):
    if VERBOSE:
        print(txt)
    out.write(txt+'\n')
    out.flush()

## Read in samples and put them in a list:
samples = []
for l in open('samples.txt'):
    if len(l) > 1:
        samples.append(l.split('/')[-1].replace('_R1_001.fastq.gz', '').strip())


# Setup folders and paths variables:
resultsDir = '01-Cleaned'
bamFolder = '02-Mapped'
variantFolder = '03-Calls'
#gatkPath = '/opt/modules/biology/gatk/3.5/bin/GenomeAnalysisTK.jar'
rawdataDir = '/mnt/lfs2/hend6746/fox_cancer/0rawdata_test'
bwaIndex = '/mnt/lfs2/hend6746/wolves/reference/canfam31/canfam31.fa'
#gatkCall = 'java -jar /opt/modules/biology/gatk/3.5/bin/GenomeAnalysisTK.jar -R %s -T HaplotypeCaller' % bwaIndex
os.system('mkdir -p %s' % resultsDir)
os.system('mkdir -p %s' % bamFolder)
os.system('mkdir -p %s' % variantFolder)

##### Run pipeline ###
for sample in samples:
    print "Processing", sample, "....."
    # Set up files:
    logFile = jp(resultsDir, sample + '_cleaning.log')
    logCommands = open(jp(resultsDir, sample + '_commands.log'), 'w')

    # First run superdeduper
    # David said don't run the compression (16Sep21)
    cmd = ' '.join(['super_deduper -1', jp(rawdataDir, sample + '_R1_001.fastq.gz'),
                    '-2', jp(rawdataDir, sample + '_R2_001.fastq.gz'), '-p', jp(resultsDir, sample + '_sd'),
                    '>>', logFile, '2>&1'])
    log(cmd, logCommands)
    os.system(cmd)

    # Second run sickle
    #The --length-threshold was set to 200, but that seems really long; default is 20
    cmd = ' '.join(['sickle pe --length-threshold 20 --qual-threshold 25 --qual-type sanger -f', jp(resultsDir, sample + '_sd_nodup_PE1.fastq'),
                    '-r', jp(resultsDir, sample + '_sd_nodup_PE2.fastq'), 
                    '--output-pe1', jp(resultsDir, sample + '_sickle_PE1.fastq'),
                    '--output-pe2', jp(resultsDir, sample + '_sickle_PE2.fastq'),
                    '--output-single', jp(resultsDir, sample + '_sickle_SE.fastq'), '>>', logFile, '2>&1'])
    log(cmd, logCommands)
    os.system(cmd)
    

    # Third run flash2
    # the --max-overlap was set to 600, but that seems really long; default is 65, try 400
    # --max-overlap 400 --min-overlap 15 --max-mismatch-density .10 --min-overlap-outie 35 --percent-cutoff 25
    cmd = ' '.join(['flash2 --max-overlap 150 --allow-outies --threads 7 -o', sample + '_flash',
                    '-d', resultsDir, jp(resultsDir, sample + '_sickle_PE1.fastq'), jp(resultsDir, sample + '_sickle_PE2.fastq'),
                     '>>', logFile, '2>&1'])
    log(cmd, logCommands)
    os.system(cmd)


    # Combine SE files:
    cmd = ' '.join(['cat', jp(resultsDir, sample + '_sickle_SE.fastq'), jp(resultsDir, sample + '_flash.extendedFrags.fastq'),
                    '>', jp(resultsDir, sample + "_cleaned_SE.fastq")])
    log(cmd, logCommands)
    os.system(cmd)

    # Rename PE and SE files to something nicer:
    cmd = ' '.join(['mv', jp(resultsDir, sample + "_flash.notCombined_1.fastq"), jp(resultsDir, sample + "_cleaned_PE1.fastq")])
    log(cmd, logCommands)
    os.system(cmd)
    cmd = ' '.join(['mv', jp(resultsDir, sample + "_flash.notCombined_2.fastq"), jp(resultsDir, sample + "_cleaned_PE2.fastq")])
    log(cmd, logCommands)
    os.system(cmd)

    # Clean up intermediary files:
    cmd = ' '.join(['rm', jp(resultsDir, sample + "_sd*"), jp(resultsDir, sample + "_sickle*"), jp(resultsDir, sample + "_flash.extendedFrags.fastq")])
    log(cmd, logCommands)
    os.system(cmd)

    # Compress cleaned files:
    cmd = ' '.join(['gzip', jp(resultsDir, '*.fastq')])
    log(cmd, logCommands)
    os.system(cmd)

    # Run BWA to map samples, combine sam files, sort
    # -t number of threads -R read group header 
    logFile = jp(bamFolder, sample + '_mapping.log')
    cmd = ' '.join(["bwa mem -t 4 -R '@RG\tID:bwa\tSM:" + sample + "\tPL:ILLUMINA'",
                    bwaIndex, jp(resultsDir, sample + "_cleaned_PE1.fastq.gz"),
                    jp(resultsDir, sample + "_cleaned_PE2.fastq.gz"), ">", jp(bamFolder, sample + "_PE.sam"),
                    "2>", logFile])
    log(cmd, logCommands)
    os.system(cmd)
    cmd = ' '.join(["bwa mem -t 4 -R '@RG\tID:bwa\tSM:" + sample + "\tPL:ILLUMINA'",
                    bwaIndex, jp(resultsDir, sample + "_cleaned_SE.fastq.gz"), ">>", jp(bamFolder, sample + "_SE.sam"),
                    "2>>", logFile])
    log(cmd, logCommands)
    os.system(cmd)

    #merge and sort
    cmd = ' '.join(['cat', jp(bamFolder, sample + "_PE.sam"), '>', jp(bamFolder, sample + ".sam")])
    log(cmd, logCommands)
    os.system(cmd)
    cmd = ' '.join(['samtools view', jp(bamFolder, sample + "_SE.sam"), '>>', jp(bamFolder, sample + ".sam")])
    log(cmd, logCommands)
    os.system(cmd)
    cmd = ' '.join(['samtools view -bS', jp(bamFolder, sample + ".sam"), '| samtools sort - -o', jp(bamFolder, sample) + ".bam"])
    log(cmd, logCommands)
    os.system(cmd)

    #Index:
    cmd = ' '.join(['samtools index', jp(bamFolder, sample) + ".bam"])
    log(cmd, logCommands)
    os.system(cmd)

    # Clean up sam files:
    cmd = ' '.join(['rm', jp(bamFolder, "*.sam")])
    os.system(cmd)
    logCommands.close()

    
'''

