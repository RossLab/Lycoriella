#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o bowtie2.$JOB_ID.log
set -e
SCRATCH=/scratch/$USER/$JOB_ID/bowtie2
DATADIR=`pwd`
mkdir -p $SCRATCH
cd $SCRATCH
# sync files in
rsync -av /data/ross/flies/raw/Lycoriella_sp/*.fq.gz "$PWD"
rsync -av /data/ross/flies/analyses/john_urban_bcop_v2_assembly/Bcop_v2-chromosomes.fasta "$PWD"

echo "building bowtie index"
bowtie2-build Bcop_v2-chromosomes.fasta Bcop_v2-chromosomes.idx

echo "starting alignment"
bowtie2 --very-sensitive-local --threads 20 \
        -x Bcop_v2-chromosomes.idx \
        -1 BCM_EDSW220012746-1a_HL5LLDSX3_L4_1.trimmed.fq.gz \
        -2 BCM_EDSW220012746-1a_HL5LLDSX3_L4_2.trimmed.fq.gz \
        | samtools view -bS - > Lyco_vs_bcopv2.bam && samtools sort Lyco_vs_bcopv2.bam \
        > Lyco_vs_bcopv2.sorted.bam

# sync files back
rsync -av *.sorted.bam /data/ross/flies/analyses/Lycoriella/results/mapping_vs_bcop/
rm *.gz *.bam *.fasta
rm -rf /scratch/$USER/$JOB_ID
