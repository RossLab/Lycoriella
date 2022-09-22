## Genome assembly

There are multiple short-read assemblers out there. Previously I have used SPAdes (https://github.com/ablab/spades). Some others that might be worth trying are:
- SOAPdenovo2 (https://github.com/aquaskyline/SOAPdenovo2)
- velvet (https://github.com/dzerbino/velvet)
- SGA (https://github.com/jts/sga)
- ALLPATHS-LG (https://github.com/danforthcenter/bioinformatics/blob/master/docs/allpaths.md)

I tried running SPAdes on the Lycoriella data but it kept crashing due to memory errors, so I tried on a subset (25%) of the reads:

```
zcat BCM_EDSW220012746-1a_HL5LLDSX3_L4_1.fq.gz | seqtk sample -s100 BCM_EDSW220012746-1a_HL5LLDSX3_L4_1.fq.gz 0.25 | gzip > BCM_subset_0.25_1.fq.gz
zcat BCM_EDSW220012746-1a_HL5LLDSX3_L4_2.fq.gz | seqtk sample -s100 BCM_EDSW220012746-1a_HL5LLDSX3_L4_2.fq.gz 0.25 | gzip > BCM_subset_0.25_2.fq.gz
rm BCM_EDSW220012746*

spades.py -t 16 --isolate -1 BCM_subset_0.25_1.fq.gz -2 BCM_subset_0.25_2.fq.gz -o lycoriella_spades/
```

This worked. You can use the program assembly-stats (https://github.com/sanger-pathogens/assembly-stats) to do a quick assessment of the genome:

```
stats for scaffolds.fasta
sum = 356096593, n = 690156, ave = 515.97, largest = 1201231
N50 = 1439, n = 33410
N60 = 837, n = 66628
N70 = 497, n = 122350
N80 = 274, n = 219160
N90 = 155, n = 392389
N100 = 78, n = 690156
N_count = 525873
Gaps = 6024
```

The things we care about here are:
- sum - number of bases or genome size
- n - number of contigs/scaffolds
- largest - size of the largest scaffold
- N50 - a bit like the median scaffold length (line up all scaffolds from shortest to longest, then N50 is the length of the scaffold 'halfway' along the genome). The larger the N50, the more contiguous the assembly.

Obviously 690k scaffolds is a crazy amount, but we expect this from a short read assembly. A lot of this will be junk, so we can filter some out using some tools seqkit (https://github.com/shenwei356/seqkit) and csvtk (https://github.com/shenwei356/csvtk). E.g. to to remove contigs with less than 5x coverage and any smaller than 500 bp (I found this command on biostars):

```
seqkit fx2tab scaffolds.fasta | csvtk mutate -H -t -f 1 -p "cov_(.+)" | csvtk mutate -H -t -f 1 -p "length_([0-9]+)" | awk -F "\t" '$4>=5 && $5>=500' | seqkit tab2fx > scaffolds.filt.fasta
```

This has got rid of a whopping ~600k scaffolds and 130Mb of sequence, but that's fine - it was probably mostly junk. And now the assembly length roughly matches up with that predicted by GenomeScope2, which is encouraging.

```
stats for scaffolds.filt.fasta
sum = 223744941, n = 95054, ave = 2353.87, largest = 1201231
N50 = 5560, n = 7615
N60 = 3406, n = 12796
N70 = 1995, n = 21430
N80 = 1208, n = 36066
N90 = 789, n = 59265
N100 = 500, n = 95054
N_count = 516720
Gaps = 5223
```

Another option is to try HGA (https://github.com/aalokaily/Hierarchical-Genome-Assembly-HGA also see the paper https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-2515-7). This overcomes the memory issue by partitioning the reads, then using SPAdes (or velvet) to assemble the subsets before combining them into a consensus sequence. It could be worth giving this a try. Unfortunately there's no conda recipe for the software, so it will need installing manually. It requires the assembly program (SPAdes) and python 2.7.6 to work. After that, it should run like this... Not exactly sure what the reads used for reassembly step are, I'll need to read up on this tool

```
python HGA.py -spades /path/to/spades -PA SPAdes -P12 reads_for_partition_step -R12 reads_for_reassembly_step -ins insert_size -P n_partitions -Rkmer kmer_size -t n_threads -out output_dir
```

## Scaffolding

If there's any significant homology between this genome and the B. coprophila genome, we could use the B. coprophila genome for scaffolding. It could be worth a try because we don't really have any other options (e.g. long reads, HiC, or a closer species) for scaffolding. RagTag (https://github.com/malonge/RagTag) is a good reference-based scaffolder.

