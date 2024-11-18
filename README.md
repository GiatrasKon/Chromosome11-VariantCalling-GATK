
# GATK-Based Variant Calling and Analysis on Chromosome 11

This repository contains a complete pipeline for germline variant calling and analysis using the Genome Analysis Toolkit (GATK). The project focuses on chromosome 11, a region of the human genome associated with numerous disease-related genes, using exome data from 27 samples from the 1000 Genomes Project. The pipeline integrates robust bioinformatics workflows for data preprocessing, variant discovery, and annotation to explore the genetic variations and their potential biological impacts.

This analysis was performed as part of the final project for the "Introduction to Bioinformatics" graduate course of the MSc Data Science & Information Technologies Master's programme (Bioinformatics - Biomedical Data Science Specialization) of the Department of Informatics and Telecommunications department of the National and Kapodistrian University of Athens (NKUA), under the supervision of professors Martin Reczko and Alexandros Dimopoulos, in the academic year 2023-2024.

### Workflow

The analysis pipeline is structured into three main phases:

1. Pre-Processing:
    - Data retrieval of chromosome 11 exome alignments.
    - BAM file quality control and alignment to the GRCh38 reference genome.
    - Removal of duplicates and base quality score recalibration (BQSR).
2. Variant Discovery:
    - Identification of single nucleotide polymorphisms (SNPs) and INDELs using GATK's HaplotypeCaller.
    - Joint genotyping of all samples to produce a cohort VCF file.
3. Callset Refinement and Annotation:
    - Variant Quality Score Recalibration (VQSR) to filter high-confidence variants.
    - Annotation using snpEff to classify variants based on genomic impact and functional categories.

### Key Results

- A total of 345,245 variants were identified, including SNPs (89.3%) and INDELs (10.7%).
- Annotation revealed most variants were located in non-coding regions (introns and intergenic regions) with potential impacts on regulatory elements.
- SNP-specific metrics, such as a Ts/Tv ratio of 2.1, align with expected genome-wide patterns.

### Cloning the Repository

```sh
git clone https://github.com/GiatrasKon/Chromosome11-VariantCalling-GATK.git
```

### Package Dependencies

Ensure you have the following packages installed:

- GATK
- SAMtools
- BWA
- bcftools
- snpEff

For package management, use `conda`:

```sh
conda create -n gatk_pipeline -c bioconda gatk samtools bwa bcftools snpeff
conda activate gatk_pipeline
```

### Documentation

Refer to the `documents` directory for the project report.

---