## Bioinformatic workflow for identification, annotation, and quantification of metagenome assembled genomes

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
```

# Co-assembly
```
megahit -1 lib1_qc_paired_fwd_read.fq.gz,lib2_qc_paired_fwd_read.fq.gz,... \
-2 lib1_qc_paired_rvs_read.fq.gz,lib2_qc_paired_rvs_read.fq.gz,... \
-o coassembly_contigs \
--out-prefix coassembly_id \
--min-contig-len 2000 \
--presets meta-large
```

# Read recruitment against assemblies
```
# Reference DB build
bowtie2-build coassembly_contigs.fa coassembly_db

# Mapping
bowtie2 -x coassembly_db \
-1 qc_paired_fwd_read.fq.gz \
-2 qc_paired_rvs_read.fq.gz \
-S coassembly_alignment.sam \
--sensitive

# Compression
samtools view -F 4 -bS coassembly_alignment.sam | samtools sort > coassembly_alignment.sI.bam
samtools index coassembly_alignment.sI.bam
```

# MAG binning
```
# Contig coverage depth quantification with MetaBAT 2
jgi_summarize_bam_contig_depths --noIntraDepthVariance --outputDepth cov_depth.txt coassembly_alignments/*.bam

# Binning with VAMB
vamb --outdir ./vamb --fasta coassembly_contigs.fa --jgi ./cov_depth.txt --minfasta 50000
```

# Dereplication
```
dRep dereplicate ./MAG_dRep \
-g ./vamb/bins/*.fna \
--S_algorithm ANImf \
-sa 0.99 \
-nc 0.1 \
-l 50000 \
-comp 50 \
-con 25 \
--checkM_method lineage_wf \
--clusterAlg single
```

# Read recruitment against MAGs
```
# Reference DB build
bowtie2-build ../drep.mags.fa mag_db

# Mapping
bowtie2 -x mag_db \
-1 qc_paired_fwd_read.fq.gz \
-2 qc_paired_rvs_read.fq.gz \
-S mag_alignment.sam \
--sensitive

# Compression
samtools view -F 4 -bS mag_alignment.sam | samtools sort > mag_alignment.sI.bam
```

# Rarefaction
```
samtools view -b -s rarefaction_factor -@ 10 mag_alignment.sI.bam | samtools sort > rare_mag_alignment.sI.bam
samtools index rare_mag_alignment.sI.bam
```

# Abundance table construction
```
coverm contig -m trimmed_mean --min-covered-fraction 0.25 -b rare_mag_alignments/*.bam > mag.tmean.tsv
```

# Taxonomic classification of MAGs
```
gtdbtk classify_wf --genome_dir ./MAG_dRep/dereplicated_genomes --out_dir gtdb

gtdbtk infer --msa_file ./gtdb/align/gtdbtk.bac120.user_msa.fasta.gz --out_dir infer_out 
```