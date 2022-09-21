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
rsync -av /data/ross/flies/raw/Lycoriella_sp/*.fq.gz "$PWD"

spades.py -t 30 --isolate -1 BCM_EDSW220012746-1a_HL5LLDSX3_L4_1.trimmed.fq.gz -2 BCM_EDSW220012746-1a_HL5LLDSX3_L4_2.trimmed.fq.gz -o lycoriella_spades/

# sync files back
rsync -av lycoriella_spades/ /data/ross/flies/analyses/Lycoriella_sp/assembly
rm *
rm -rf /scratch/$USER/$JOB_ID
