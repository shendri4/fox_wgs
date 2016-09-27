#!/usr/bin/env python	
###	Usage ./fastqc.py

### Make sure to load:
### module load fastqc	

from os.path import join as jp
import os
import sys

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
rawdataDir = '/mnt/lfs2/hend6746/fox_cancer/0rawdata_test'
cleandataDir = '01-Cleaned'
bamFolder = '02-Mapped'
resultsDir = '/mnt/lfs2/hend6746/fox_cancer/0rawdata_test/fastqc_test'
os.system('mkdir -p %s' % resultsDir)


##### Run pipeline ###
for sample in samples:
    print "Processing", sample, "....."
    # Set up files:
    logFile = jp(resultsDir, sample + '_cleaning.log')
    logCommands = open(jp(resultsDir, sample + '_commands.log'), 'w')

    #raw
    cmd = ' '.join(['fastqc', '--outdir', resultsDir, '--format fastq', jp(rawdataDir, sample + '_R1_001.fastq.gz'), '>>', logFile, '2>&1']) 
    os.system(cmd)
    logCommands.close()
    
    cmd = ' '.join(['fastqc', '--outdir', resultsDir, '--format fastq', jp(rawdataDir, sample + '_R2_001.fastq.gz'), '>>', logFile, '2>&1'])
    os.system(cmd)
    logCommands.close()
'''    
    #cleaned
    cmd = ' '.join(['fastqc', '--outdir', resultsDir, '--format', fastq, jp(cleandataDir, sample + '_cleaned_PE1.fastq.gz')]) 
    cmd = ' '.join(['fastqc', '--outdir', resultsDir, '--format', fastq, jp(cleandataDir, sample + '_cleaned_PE2.fastq.gz')])
    cmd = ' '.join(['fastqc', '--outdir', resultsDir, '--format', fastq, jp(cleandataDir, sample + '_cleaned_SE.fastq.gz')]) 

    #bam
    cmd = ' '.join(['fastqc', '--outdir', resultsDir, '--format', bam, jp(bamFolder, sample + '.bam')]) 
'''