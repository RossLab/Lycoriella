#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o bwa.$JOB_ID.log
set -e
SCRATCH=/scratch/$USER/$JOB_ID/bwa
DATADIR=`pwd`
mkdir -p $SCRATCH
cd $SCRATCH
# sync files in
rsync -av /data/ross/flies/raw/Lycoriella_sp/*.fq.gz "$PWD"
rsync -av /data/ross/flies/analyses/john_urban_bcop_v2_assembly/Bcop_v2-chromosomes.fasta "$PWD"

bwa index Bcop_v2-chromosomes.fasta
bwa mem Bcop_v2-chromosomes.fasta BCM_EDSW220012746-1a_HL5LLDSX3_L4_1.trimmed.fq.gz BCM_EDSW220012746-1a_HL5LLDSX3_L4_2.trimmed.fq.gz | samtools view -b - > Lyco_vs_bcopv2.bam && samtools sort -@ 8 -o Lyco_vs_bcopv2.sorted.bam Lyco_vs_bcopv2.bam

# sync files back
rsync -av *.sorted.bam /data/ross/flies/analyses/Lycoriella_sp
rm *
rm -rf /scratch/$USER/$JOB_ID
