## Bioinformatic workflow for identification, annotation, and quantification of viral operational taxonomic units

# Quality filtering
```
#Trimmomatic
trimmomatic PE -threads 8 -phred33 \
  input_fwd_read.fq.gz input_rvs_read.fq.gz \
  trimmed_paired_fwd_read.fq.gz trimmed_unpaired_fwd_read.fq.gz \
  trimmed_paired_rvs_read.fq.gz trimmed_unpaired_rvs_read.fq.gz \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
  SLIDINGWINDOW:4:30 MINLEN:50
  
#BBDuk on paired reads
bbduk.sh \
  in1=trimmed_paired_fwd_read.fq.gz \
  in2=trimmed_paired_rvs_read.fq.gz \
  out1=qc_paired_fwd_read.fq.gz \
  out2=qc_paired_rvs_read.fq.gz \
  ref=phix174_ill.ref.fa \
  k=31 \
  hdist=1 \
  stats=stats.txt 
  
#BBDuk on unpaired reads
bbduk.sh \
    in=trimmed_unpaired_fwd_read.fq.gz \
    out=qc_unpaired_fwd_read.fq.gz \
    ref=phix174_ill.ref.fa \
    k=31 \
    hdist=1 \
    stats=stats.txt
```

# Assembly
```
megahit -1 qc_paired_fwd_read.fq.gz \
-2 qc_paired_rvs_read.fq.gz \
-r qc_unpaired_fwd_read.fq.gz,qc_unpaired_rvs_read.fq.gz \
-o contigs \
--out-prefix sample_id \
--min-contig-len 10000 \
--presets meta-large
```

# Viral classification
```
VIBRANT_run.py -i megahit.contigs.fa \
-folder vibrant \
-f nucl \
-virome #remove virome tag when processing total metagenomes
```

# Dereplication
```
dRep dereplicate dRep \
-g ./viral_contigs/*.fa \
--S_algorithm ANImf \
-sa 0.95 \
-nc 0.85 \
-l 10000 \
-N50W 0 \
-sizeW 1 \
--ignoreGenomeQuality \
--clusterAlg single 
```

# Read recruitment to dereplicated vOTUs
```
# Reference DB build
bowtie2-build ../drep.votus.fa votu_db

# Mapping
bowtie2 -x votu_db \
-1 qc_paired_fwd_read.fq.gz \
-2 qc_paired_rvs_read.fq.gz \
-S votu_alignment.sam \
--sensitive

# Compression
samtools view -F 4 -bS votu_alignment.sam | samtools sort > votu_alignment.sI.bam
```

# Rarefaction
```
samtools view -b -s rarefaction_factor -@ 10 votu_alignment.sI.bam | samtools sort > rare_votu_alignment.sI.bam
samtools index rare_votu_alignment.sI.bam
```

# Abundance table construction
```
coverm contig -m trimmed_mean --min-covered-fraction 0.75 -b rare_votu_alignments/*.bam > votu.tmean.tsv
```

# Gene prediction
```
prodigal -i drep.votus.fa -a votus.prodigal.faa -p meta
```

# VCONTACT
```
vcontact2 --raw-proteins votus.prodigal.faa \
--rel-mode 'Diamond' \
--db 'ProkaryoticViralRefSeq85-Merged' \
--proteins-fp gene2genome.csv \
--pcs-mode MCL \
--vcs-mode ClusterONE \
--c1-bin cluster_one-1.0.jar \
--output-dir vcontact_out
```

# Host prediction
```
iphop predict --fa_file drep.votus.fa \
--out_dir ./iphop_results \
--db_dir iphop_db
```

