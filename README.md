****M. leprae WGS Analysis Pipeline****

A Bash-script pipeline for whole-genome sequencing (WGS) analysis of Mycobacterium leprae, covering quality control, read processing, alignment, variant calling, filtering, and annotation.

**Workflow**

SRA → FASTQ FILE → QC → Trimming → BWA Alignment
→ BAM Processing → Variant Calling → Filtering
→ VCF Normalization → SnpEff Annotation

**Tools**

SRA Toolkit,
FastQC,
fastp,
BWA,
SAMtools,
Picard,
BCFtools,
SnpEff,
SnpSift,
Qualimap,
MultiQC,

**Repository Structure**

M.-leprae-pipeline/

├── README.md

├── pipeline.sh

└── wgs_env.yml

**Installation**

git clone https://github.com/tusharbioinfo/M.-leprae-pipeline.git

cd M.-leprae-pipeline

conda env create -f wgs_env.yml

conda activate wgs

chmod +x pipeline.sh

Input-

The pipeline requires a CSV file containing SRA accession IDs and sequencing type.

sra_id,type

SRRXXXXXXX,PAIRED
SRRXXXXXXX,PAIRED
SRRXXXXXXX,SINGLE

Run

./pipeline.sh samples.csv

Analysis-
The pipeline performs:

SRA data download and FASTQ conversion
Read quality control
Adapter and quality trimming
Reference genome alignment
BAM processing and sorting
Duplicate removal
Variant calling and filtering
VCF normalization
Variant annotation
Quality assessment and reporting

Reference Genome-
Mycobacterium leprae reference genome:
GCF_000195855.1

Applications-
Genomic diversity analysis
SNP identification
Variant annotation
Comparative genomics
Antimicrobial-resistance-associated variant analysis
Downstream phylogenetic analysis
Author

Tushar Sain
