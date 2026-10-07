M. leprae WGS Analysis Pipeline
A Bash script pipeline for whole-genome sequencing (WGS) analysis of Mycobacterium leprae, covering quality control, read processing, alignment, variant calling, filtering, and annotation.
Workflow
SRA → FASTQ FILE → QC → Trimming → BWA Alignment
→ BAM Processing → Variant Calling → Filtering
→ VCF Normalization → SnpEff Annotation
Tools
•	SRA Toolkit
•	FastQC
•	fastp
•	BWA
•	SAMtools
•	Picard
•	BCFtools
•	SnpEff
•	SnpSift
•	Qualimap
•	MultiQC
Repository Structure
M.-leprae-pipeline/
├── README.md
├── pipeline.sh
└── wgs_env.yml
Installation
git clone https://github.com/tusharbioinfo/M.-leprae-pipeline.git
cd M.-leprae-pipeline

conda env create -f wgs_env.yml
conda activate wgs

chmod +x pipeline.sh
Input
The pipeline requires a CSV(text) file containing SRA accession IDs and sequencing type.
sra_id,type
SRRXXXXXXX,PAIRED
SRRXXXXXXX,PAIRED
SRRXXXXXXX,SINGLE
Run
./pipeline.sh samples.csv
Analysis
The pipeline performs:
1.	SRA data download and FASTQ conversion
2.	Read quality control
3.	Adapter and quality trimming
4.	Reference genome alignment
5.	BAM processing and sorting
6.	Duplicate removal
7.	Variant calling and filtering
8.	VCF normalization
9.	Variant annotation
10.	Quality assessment and reporting
Reference Genome
Mycobacterium leprae reference genome:
GCF_000195855.1
Applications
•	Genomic diversity analysis
•	SNP identification
•	Variant annotation
•	Comparative genomics
•	Antimicrobial-resistance-associated variant analysis
•	Phylogenetic analysis
Author
Tushar Sain
Bioinformatics | Genomics | NGS | Computational Biology
