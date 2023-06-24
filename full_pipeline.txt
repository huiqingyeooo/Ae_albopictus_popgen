################################
# Full pipeline (24 June 2023)
################################

# Quality check of raw reads
fastqc *.fastq.gz --outdir=./1_sequences

# Filter low quality sequences
for f in $(cat ./list.txt)
do
cutadapt -q 30 -o ${f}_q30.fastq.gz ./${f}.fastq.gz
done

# Demultiplexing sequence reads
process_radtags -1 ./pool_1.fq.gz -2 ./pool_2.fq.gz -o ./2_process_radtags/pool --renz_1 ecoRI --renz_2 mspI -b ./barcodes_pool.txt --inline_index -c -q -r

# Repeat mask reference genome
bedtools maskfasta -fi genome.fasta -bed repeat.gff -fo masked_genome.fasta

# Index reference genome
bwa index ./ GCA_001444175.1_A.albopictus_v1.0_genomic_RepeatMasked.fasta

# Rename headers in sequence file
for f in $(cat ./list.txt)
do
./bbrename.sh in1=./${f}.1.fq.gz in2=./${f}.2.fq.gz \
out1=./output/${f}.1_renamed.fq.gz out2=./output/${f}.2_renamed.fq.gz \
prefix=read
done

# Alignment to reference genome
for f in $(cat ./list.txt)
do
bwa mem -t 8 -M ./GCA_001444175.1_A.albopictus_v1.0_genomic_RepeatMasked.fasta \
./output/${f}.1_renamed.fq.gz ./output/${f}.2_renamed.fq.gz > ${f}.sam
done

# Convert aligned file from .sam to .bam
for f in $(cat ./list.txt)
do
samtools view -bShq 20 ./${f}.sam >  ./${f}.bam
done

# Alignment statistics
for f in $(cat ./list.txt)
do
samtools flagstat -@15 ./bam_files_q20/${f}.bam > ./_flagstat/${f}.txt
done

# Sort .bam files
for f in $(cat ./list.txt)
do
samtools sort -m 2G -@ 15 -o ./bam_files_q20_sorted/${f}.bam ./bam_files_q20/${f}.bam
done

# SNP calling with Stacks
ref_map.pl -T 15 --samples ./bam_files_q20_sorted -o ./6_refmap --popmap ./6_refmap.txt

populations -P ./6_refmap -O ./bam_files_q20_sorted -M ./pops.txt -t 15 -r 0.75 --min-maf 0.05 --write-single-snp --fstats --hwe --structure --genepop --plink --vcf --fasta -k –verbose

# Use plink to filter SNP set
  ## indep pairwise-linkage
  plink –file ./populations.plink --allow-no-sex --allow-extra-chr --indep-pairwise 25 10 0.9

  ## Extract SNP data
  plink --file ./populations.plink --extract plink.prune.in --allow-extra-chr --make-bed –out albo_prune

  ## Format output in vcf format, remove high missing samples
  plink --bfile albo_prune --allow-extra-chr --recode vcf --out albo_plink

  ## Format output in structure format
  plink --bfile albo _prune --allow-extra-chr --recode structure --out albo_plink

  ## Check missing data for each individual
  plink --bfile albo_prune --allow-extra-chr --missing --allow-no-sex –out albo_postprune

# Calculate basic popgen statistics (e.g., heterozygosity, inbreeding coefficient)
Please refer to heterozygosity.R

# Identify close kins and remove one individual from close kin pairs
Please refer to identity_by_descent.R

# Principal Component Analysis with SNPRelate
Please refer to PCA.R

# Test for isolation-by-distance with Mantel test
Please refer to partial_mantel_test.R

# Admixture analysis
  ## Make first column a bunch of zeros
  awk '{$1=0;print $0}' ./albo_prune.bim > ./albo_prune_format.bim
  rm albo_prune.bim
  mv albo_prune_format.bim albo_prune.bim

  ## Run admixture
  admixture ./albo_prune.bed 2

  ## To choose the correct value for k, use ADMIXTURE cross validation procedure
  for K in 1 2 3 4 5
  do
  admixture --cv albo_prune.bed $K | tee log${K}.out
  done

  ## To view CV errors. Lowest error = best K value
  grep -h CV log*.out > CV_error.txt

# STRUCTURE analysis
structure_threader run -K 2 -R 4 -i ./albo -o ./albo_structure -t 8 -st ./structure --params ./mainparams

# Demographic history analysis
  ## SNP calling with angsd
  ./angsd -bam ./bam.filelist -out ./sfs -anc ./GCA_001444175.1_A.albopictus_v1.0_genomic_RepeatMasked.fasta -GL 1 -doCounts 1 -dosaf 1 -doDepth 1 -nThreads 8 -    dumpCounts 1 -minMapQ 20 -minQ 20

  ## Obtain site frequency spectrum
  ./angsd/misc/realSFS -P 8 ./sfs.saf.idx -fold 1 > ./folded.sfs

  ## Stairway plot2
  java -cp stairway_plot_es Stairbuilder albo.blueprint
  bash albo.blueprint.sh

  ## Refer to plot.stairwayplot2.R for code to generate  graph

# Landscape genetic anaylsis with ResDisMapper
Please refer to ResDisMapper.R

# Landscape genetic optimization with ResistanceGA
Please refer to ResistanceGA.R

# Wolbachia infection analysis
Please refer to Wolbachia_analysis.R
