#=====================================PRE-PROCESSING========================================

#=====================================Data Retrieval========================================
# Static parts of the URL
link_prefix="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/"
link_suffix_1=".chrom11.ILLUMINA.bwa.GBR.exome.20120522.bam"
link_suffix_2=".chrom11.ILLUMINA.bwa.GBR.exome.20121211.bam"
link_middle="/exome_alignment/"

# Static sample prefix
sample_prefix="HG00"

# wget query for all the possible link combinations
for i in $(seq -w 100 140); do
  wget ${link_prefix}${sample_prefix}${i}${link_middle}${sample_prefix}${i}${link_suffix_1}
  wget ${link_prefix}${sample_prefix}${i}${link_middle}${sample_prefix}${i}${link_suffix_2}
done

#===============BAM File Quality Control and Reference Genome Retrieval=====================
for file in *.bam; do samtools flagstat $file | grep "QC";done # checking the quality of the BAM files

for file in *.bam; do echo $file; samtools view -H $file | grep '^@SQ'| grep 'UR';done # Checking the reference assembly identifier

samtools faidx hg38.fa chr11 > chr11.fa # Extracting the chromosome 11 from the reference genome
samtools faidx chr11.fa # Indexing the chromosome 11 reference genome

gatk CreateSequenceDictionary R=my_sequence.fa O=my_sequence.dict # Creating a sequence dictionary for the reference genome

#================Realignment to Reference Genome GRCh38, Chromosome 11======================
# Loop converting BAM files to fastq format and shuffling reads
echo "Beginning conversion from BAM to fastq"
for file in *.bam; do
  echo "Starting conversion of $file"
  # Using samtools collate to shuffle and convert BAM to fastq, then compress with gzip
  samtools collate -uOn 128 "$file" "tmp_${file}" | samtools fastq | gzip >  "interleaved_${file%.bam}.fq.gz" 
  echo "Conversion of $file complete"
done

bwa index chr11.fa # Indexing the reference genome

echo "Beginning alignment"
for file in H*.bam; do
  echo "Starting alignment of $file"
  # Extracting the sample name by removing the expected patterns from the BAM filename
  samplename=$(echo $file | sed -E 's/.chrom11.ILLUMINA.bwa.GBR.exome.(20120522|20121211).bam//')

  # Constructing the read group information string with the -R option
  read_group="@RG\tID:${samplename}\tLB:${samplename}\tSM:${samplename}\tPL:ILLUMINA\tCN:unknown"

  # Looking for the FASTQ files with both possible dates
  fq1=$(ls interleaved_${samplename}.chrom11.ILLUMINA.bwa.GBR.exome.20120522.fq.gz 2> /dev/null)
  fq2=$(ls interleaved_${samplename}.chrom11.ILLUMINA.bwa.GBR.exome.20121211.fq.gz 2> /dev/null)

  # Selecting the available FASTQ file
  if [[ -n $fq1 ]]; then
    fastq=$fq1
  elif [[ -n $fq2 ]]; then
    fastq=$fq2
  else
    echo "No FASTQ file found for $samplename"
    continue
  fi

  echo "Running BWA MEM for $samplename with FASTQ file $fastq"
  bwa mem -t 7 -M -R "$read_group" chr11.fa $fastq | samtools view -bS - | samtools sort -o "aligned${samplename}.bam"
  echo "Alignment of $file complete"
done

# This script handles the collection of total and mapped reads before and after the alignment
touch old_mapped_stats.csv
echo "sample,old bam mapped,old bam reads" >> old_mapped_stats.csv
for file in H*.bam;
do
  sample=${file:0:7}
  previous_unmapped_reads=$(samtools view -F 4 -c $file)
  previous_total_reads=$(samtools view -c $file)
  echo "$sample,$previous_unmapped_reads,$previous_total_reads" >> old_mapped_stats.csv
done

touch new_mapped_stats.csv
echo "new bam mapped,new bam reads" >> new_mapped_stats.csv
for file in aligned*.bam;
do
  reference_unmapped_reads=$(samtools view -F 4 -c $file)
  reference_aligned_reads=$(samtools view -c $file)
  echo "$reference_unmapped_reads,$reference_aligned_reads" >> new_mapped_stats.csv
done
# Pasting the files together
paste -d "," old_mapped_stats.csv new_mapped_stats.csv > mapped_stats.csv

#=====================================Mark Duplicates=======================================
# Locating and tag duplicate reads in the aligned bam files, as well as create an index bam file
for file in aligned*.bam; do
  echo "Beginning duplicate marking of $file"
  gatk MarkDuplicates \
    I="$file" \
    O="marked_dup_${file}" \
    M="marked_dup_metrics_${file%.bam}.txt" \
    CREATE_INDEX=true
  echo "Duplicate marking of $file complete"
done

# Creating a file to store duplication percentage
echo "Sample,Duplication_Percentage" > duplicate_percentage.csv
for file in marked_dup_metrics_*.txt; do
  sample_name=$(basename "$file" .txt | sed 's/marked_dup_metrics_//')
  duplication_percentage=$(awk 'NR>7 {print $9}' "$file")
  echo "$sample_name,$duplication_percentage" >> duplicate_percentage.csv
done

#=========================Base Quality Score Recalibration (BQSR)===========================
# First pass - Creating recalibration tables for each marked duplicate BAM file
for file in marked_dup_aligned*.bam; do
  echo "Creating recalibration table corresponding to $file"
  recal_table="recal_${file%.bam}.table"
  gatk BaseRecalibrator \
    -R chr11.fa \
    -I "$file" \
    --known-sites Homo_sapiens_assembly38.dbsnp138.vcf \
    -O "$recal_table"
  echo "Recalibration table created for $file"
done

# Data CSV retrieval for each recalibration table
for table in recal_*table; do
  echo "Analyzing covariates for $table"
  csv_file="../statistics/${table%.table}.csv"
  gatk AnalyzeCovariates \
    -bqsr "$table" \
    -plots "$csv_file"
  echo "Covariate analysis complete for $table, CSV file created at $csv_file"
done

# Second pass - Applying recalibration to each marked duplicate BAM file
for file in marked_dup_aligned*.bam; do
  echo "Applying recalibration to $file"
  recal_table="recal_${file%.bam}.table"
  output_bam="recalibrated_${file}"
  gatk ApplyBQSR \
    -R chr11.fa \
    -I "$file" \
    --bqsr-recal-file "$recal_table" \
    -O "$output_bam"
  echo "Recalibration applied to $file"
done

#===================================VARIANT DISCOVERY=======================================

#====================================HaplotypeCaller========================================
# Creation of GVCF files using GATK HaplotypeCaller
for file in recalibrated_*.bam; do
  echo "Starting creation of GVCF for $file"
  output_vcf="../vcf/${file##*/}"
  output_vcf="${output_vcf%.bam}.g.vcf"
  gatk HaplotypeCaller \
    -R chr11.fa \
    -I "$file" \
    -O "$output_vcf" \
    --native-pair-hmm-threads 7 \
    -ERC GVCF \
    -L chr11
  echo "Created GVCF for $file"
done

# Collection of summary statistics from the GVCF files and appending to a single CSV file
stats_file="../statistics/HC_stats.csv"

echo "Sample,number of samples,number of records,number of no-ALTs,number of SNPs,number of MNPs,number of indels,number of others,number of multiallelic sites,number of multiallelic SNP sites" > "$stats_file"

for file in ../vcf/*.g.vcf; do
  sample_name=$(basename "$file" .g.vcf)
  stats=$(bcftools stats "$file" | grep "^SN" | head -n 9 | awk -F':' '{printf("%s,", $2)}' | sed 's/,$//')
  echo "$sample_name,$stats" >> "$stats_file"
done

#====================================Consolidate GVCFs======================================
# Creating the sample map for GenomicsDBImport
cohort_sample_map="../vcf/cohort.sample_map"
touch "$cohort_sample_map"

# Looping through the g.vcf files in the 'vcf' directory to create the sample map
for file in ../vcf/*.g.vcf; do
  sample=$(basename "$file" .g.vcf)
  echo -e "$sample\t$file" >> "$cohort_sample_map"
done

# Executing GenomicsDBImport
gatk GenomicsDBImport \
  --genomicsdb-workspace-path ../database \
  --sample-name-map "$cohort_sample_map" \
  --reader-threads 7 \
  -L chr11
echo "GenomicsDBImport complete"

#=====================================Genotype GVCFs========================================
# GenotypeGVCFs
echo "Running GenotypeGVCFs"
gatk GenotypeGVCFs \
  -R chr11.fa \
  -V gendb://../database \
  -O ../vcf/chr11_cohort.vcf \
  -L chr11
echo "GenotypeGVCFs complete"

# Retrieving metrics for the cohort VCF
cohort_stats_file="../statistics/cohort_stats.csv"
echo "cohort_file,number of samples,number of records,number of no-ALTs,number of SNPs,number of MNPs,number of indels,number of others,number of multiallelic sites,number of multiallelic SNP sites" > "$cohort_stats_file"

echo "Processing chr11_cohort.vcf"
stats=$(bcftools stats ../vcf/chr11_cohort.vcf | grep "^SN" | head -n 9 | awk -F':' '{printf("%s,", $2)}' | sed 's/,$//')
echo "chr11_cohort.vcf,$stats" >> "$cohort_stats_file"
echo "Metrics for chr11_cohort.vcf saved to $cohort_stats_file"

#=====================================CALLSET REFINEMENT=====================================

#============================Variant Quality Score Recalibration=============================
# Beginning SNP Filtering
echo "Beginning SNP Filtering"
gatk VariantRecalibrator \
  -R chr11.fa \
  -V ../vcf/chr11_cohort.vcf \
  --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz \
  --resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.vcf.gz \
  --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 Homo_sapiens_assembly38.dbsnp138.vcf \
  -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
  -mode SNP \
  -O ../vcf/SNP.recal \
  --tranches-file ../statistics/SNP.tranches \
  --rscript-file ../statistics/SNP.plots.R

# Beginning INDEL Filtering
echo "Beginning INDEL Filtering"
gatk VariantRecalibrator \
  -R chr11.fa \
  -V ../vcf/chr11_cohort.vcf \
  --resource:mills,known=false,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --resource:indels,known=true,training=false,truth=false,prior=2.0 Homo_sapiens_assembly38.known_indels.vcf.gz \
  -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
  -mode INDEL \
  -O ../vcf/INDEL.recal \
  --tranches-file ../statistics/INDEL.tranches \
  --rscript-file ../statistics/INDEL.plots.R

# Applying SNP filters
echo "Applying SNP filters"
gatk ApplyVQSR \
  -R chr11.fa \
  -V ../vcf/chr11_cohort.vcf \
  -O ../vcf/SNP_filtered_chr11_cohort.vcf \
  --tranches-file ../statistics/SNP.tranches \
  --recal-file ../vcf/SNP.recal \
  --truth-sensitivity-filter-level 99.0 \
  -mode SNP

# Applying INDEL filters
echo "Applying INDEL filters"
gatk ApplyVQSR \
  -R chr11.fa \
  -V ../vcf/SNP_filtered_chr11_cohort.vcf \
  -O ../vcf/recal_filtered_chr11_cohort.vcf \
  --tranches-file ../statistics/INDEL.tranches \
  --recal-file ../vcf/INDEL.recal \
  --truth-sensitivity-filter-level 99.0 \
  -mode INDEL

#==================================Relevant Sample Isolation=================================
# Creating a relevant_samples.txt file
cat << EOF > relevant_samples.txt
recalibrated_marked_dup_alignedHG00100
recalibrated_marked_dup_alignedHG00101
recalibrated_marked_dup_alignedHG00102
recalibrated_marked_dup_alignedHG00103
recalibrated_marked_dup_alignedHG00107
recalibrated_marked_dup_alignedHG00108
recalibrated_marked_dup_alignedHG00109
recalibrated_marked_dup_alignedHG00110
recalibrated_marked_dup_alignedHG00111
recalibrated_marked_dup_alignedHG00113
recalibrated_marked_dup_alignedHG00114
EOF

# Filtering out the relevant samples from the VCF file
echo "Filtering relevant samples from recalibrated VCF"
bcftools view --samples-file relevant_samples.txt ../vcf/recal_filtered_chr11_cohort.vcf > ../vcf/relevant_recal_filtered_chr11_cohort.vcf

#=====================================Variant Annotation======================================
# Running snpEff for variant annotation
echo "Annotating filtered VCF with snpEff"
java -jar ../../../apps/snpEff/snpEff.jar -v -stats ../statistics/annot_metrics.html -csvStats ../statistics/annot_metrics.csv hg38 relevant_recal_filtered_chr11_cohort.vcf > annotated_relevant_recal_filtered_chr11_cohort.ann.vcf