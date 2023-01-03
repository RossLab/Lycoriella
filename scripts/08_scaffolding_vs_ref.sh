#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o ragtag.$JOB_ID.log
set -e
SCRATCH=/scratch/$USER/$JOB_ID/ragtag
DATADIR=`pwd`
mkdir -p $SCRATCH
cd $SCRATCH
# sync files in
rsync -av /data/ross/flies/analyses/john_urban_bcop_v2_assembly/Bcop_v2-chromosomes.fasta "$PWD"
rsync -av /data/ross/flies/analyses/Lycoriella/results/genome_assembly/lycoriella_spades/scaffolds.fasta "$PWD"

# scaffolding...
ragtag.py scaffold -t 16 Bcop_v2-chromosomes.fasta scaffolds.fasta
