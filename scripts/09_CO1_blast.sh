#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o blastn.$JOB_ID.log
set -e
SCRATCH=/scratch/$USER/$JOB_ID/blastn
DATADIR=`pwd`
mkdir -p $SCRATCH
cd $SCRATCH

rsync -av /data/ross/flies/analyses/bradysia_sex_determination/021_phylogeny/sequences/Sciaridae_with_known_strategies_CO1s.fasta "$PWD"
rsync -av /data/ross/flies/analyses/Lycoriella/results/genome_assembly/lycoriella_spades/scaffolds.filt.fasta "$PWD"

makeblastdb -in scaffolds.filt.fasta -dbtype nucl

blastn -query Sciaridae_with_known_strategies_CO1s.fasta -db scaffolds.filt.fasta -max_target_seqs 10 -evalue 1e-10 -outfmt "6 std qlen slen" -num_threads 16 > Ling_CO1.blast

rsync -av Ling_CO1.blast /data/ross/flies/analyses/Lycoriella/results/blast

rm -rf /scratch/$USER/$JOB_ID
