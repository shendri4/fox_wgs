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

#print samples

# Setup folders and paths variables:
resultsDir = '01-Cleaned'
bamFolder = '02-Mapped'
variantFolder = '03-Calls'
gatkPath = '/opt/modules/biology/gatk/3.5/bin/GenomeAnalysisTK.jar'
rawdataDir = '/mnt/lfs2/hend6746/fox_cancer/0rawdata_test'
bwaIndex = '/mnt/lfs2/hend6746/wolves/reference/canfam31/canfam31.fa'
gatkCall = 'java -jar /opt/modules/biology/gatk/3.5/bin/GenomeAnalysisTK.jar -R %s -T HaplotypeCaller' % bwaIndex
os.system('mkdir -p %s' % resultsDir)
os.system('mkdir -p %s' % bamFolder)
os.system('mkdir -p %s' % variantFolder)

##### Run pipeline ###
for sample in samples:
    print "Processing", sample, "....."
    # Set up files:
    logCommands = open(jp(resultsDir, sample + '_commands.log'), 'w')

    #Call SNPs with GATK
    logFile = jp(variantFolder, sample + '_GATK.log')
    cmd = ' '.join([gatkCall,  ' -I ' + jp(bamFolder, sample) + ".bam", 
                    ' -o ' + jp(variantFolder, sample) + ".raw.variants.vcf", '>>', logFile, '2>&1'])
    log(cmd, logCommands)
    os.system(cmd)
    
'''
gatkCall += ' -I ' + jp(bamFolder, sample) + ".bam"
logCommands.close()
gatkCall += ' -o output.raw.variants.vcf'

print "Now call gatk with:\n"+gatkCall
with open("run_gatk.sh", 'w') as outf:
    outf.write(gatkCall + '\n')
    
##### Merge all bam files
#cmd = ' '.join(['samtools merge', jp(bamFolder, "*.bam")])
'''