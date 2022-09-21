#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o spades.$JOB_ID.log
set -e
SCRATCH=/scratch/$USER/$JOB_ID/spades
DATADIR=`pwd`
mkdir -p $SCRATCH
cd $SCRATCH
# sync files in
rsync -av /data/ross/sequencing/raw/novogene_june_2022_mealybugs_gnats/raw_data/BCM/*.fq.gz "$PWD"

zcat BCM_EDSW220012746-1a_HL5LLDSX3_L4_1.fq.gz | seqtk sample -s100 BCM_EDSW220012746-1a_HL5LLDSX3_L4_1.fq.gz 0.25 | gzip > BCM_subset_0.25_1.fq.gz
zcat BCM_EDSW220012746-1a_HL5LLDSX3_L4_2.fq.gz | seqtk sample -s100 BCM_EDSW220012746-1a_HL5LLDSX3_L4_2.fq.gz 0.25 | gzip > BCM_subset_0.25_2.fq.gz
rm BCM_EDSW220012746*

spades.py -t 30 --isolate -1 BCM_subset_0.25_1.fq.gz -2 BCM_subset_0.25_2.fq.gz -o lycoriella_spades/

# sync files back
rsync -av lycoriella_spades /data/ross/flies/analyses/Lycoriella/results/genome_assembly
rm *.fq*
rm lycoriella_spades/*
rm -rf /scratch/$USER/$JOB_ID
