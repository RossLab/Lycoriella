#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o qc_trim.$JOB_ID.log
set -e
SCRATCH=/scratch/$USER/$JOB_ID/qc_trim
DATADIR=`pwd`
mkdir -p $SCRATCH
cd $SCRATCH
# sync files in
rsync -av /data/ross/sequencing/raw/novogene_june_2022_mealybugs_gnats/raw_data/BCM/*.gz "$PWD"

# fastqc pre-trim
fastqc -t 4 *.fq.gz

# trim reads
for file in $(ls *.fq.gz)
do
	base=$(basename $file "_1.fq.gz")
  	fastp -i ${base}_1.fq.gz -I ${base}_2.fq.gz -o ${base}_1.trimmed.fq.gz -O ${base}_2.trimmed.fq.gz
done

# fastqc post-trim
fastqc -t 4 *.trimmed.fq.gz

rsync -av *.html /data/ross/flies/analyses/Lycoriella_sp/qc
rsync -av *.trimmed.fq.gz /data/ross/flies/raw/Lycoriella_sp/
rm *.gz
rm *.html
rm -rf /scratch/$USER/$JOB_ID
