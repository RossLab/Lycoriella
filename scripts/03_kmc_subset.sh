#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o kmc.$JOB_ID.log
set -e
SCRATCH=/scratch/$USER/$JOB_ID/kmc
DATADIR=`pwd`
mkdir -p $SCRATCH
cd $SCRATCH
# sync files in
rsync -av /data/ross/sequencing/raw/novogene_june_2022_mealybugs_gnats/raw_data/BCM/*.fq.gz "$PWD"

zcat BCM_EDSW220012746-1a_HL5LLDSX3_L4_1.fq.gz | seqtk sample -s100 BCM_EDSW220012746-1a_HL5LLDSX3_L4_1.fq.gz 0.1 | gzip > BCM_subset_30m_1.fq.gz
zcat BCM_EDSW220012746-1a_HL5LLDSX3_L4_2.fq.gz | seqtk sample -s100 BCM_EDSW220012746-1a_HL5LLDSX3_L4_2.fq.gz 0.1 | gzip > BCM_subset_30m_2.fq.gz
rm BCM_EDSW220012746*

for file in *_1.fq.gz
do
	base=$(basename $file "_1.fq.gz")
	mkdir ${base}_tmp/
	cat ${base}_1.fq.gz ${base}_2.fq.gz > ${base}_files.fq.gz
	kmc -k21 -t30 -m64 -ci1 -cs10000 -fq ${base}_files.fq.gz ${base}_kmer_counts ${base}_tmp/
	kmc_tools transform ${base}_kmer_counts histogram ${base}_0.1_subset_kmer_k21.hist -cx100000
done

# sync files back
rsync -av *.hist /data/ross/flies/analyses/Lycoriella_sp/kmc
rm *.fq*
rm -rf /scratch/$USER/$JOB_ID
