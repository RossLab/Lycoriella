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
rsync -av /data/ross/flies/analyses/Lycoriella_sp/mapping/*.bam "$PWD"

for file in $(ls *.bam)
do
	base=$(basename $file ".bam")
	samtools index ${base}.bam
	bedtools genomecov -ibam ${base}.bam -d > per_base_cov.txt
	Rscript /data/ross/flies/analyses/Lycoriella_sp/scripts/get_100kb_cov_windows.R
	mv 100kb_cov_win_means.txt ${base}.100kb.win.means.txt
done

rsync -av *.100kb.win.means.txt /data/ross/flies/analyses/Lycoriella_sp/mapping/
rm *
rm -rf /scratch/$USER/$JOB_ID
