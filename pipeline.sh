#!/bin/bash
set -o pipefail
# A simple pipeline script for M leprae NGS data analysis starting from data download 
# Author: Arup Ghosh, Tushar Sain
# Date created: 17-04-2024
# Last modified: 17-04-2024
# Exclude the header line

#downolad refernce genome (TN)
#bwa_index= ../references/m_leprae/bwa/leprae
#samtools faidx ./index/GCF_000195855.1_ASM19585v1_genomic.fna

#create folders for pipeline run
mkdir -p fastq/trimmed fastq/rawdata reports/fastqc_raw results/aligned_bam results/sorted_bam results/dedup_bam results/sorted_bam results/varient_vcf results/pileup results/filtered_vcf reports/fastqc_trimmed reports/samtools_sort reports/samtools_dedup index results/anotation reports/dedup_stat reports/qualimap
#downolad refernce genome (TN)
curl "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000195855.1/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED" --output ./index/genome.zip

unzip  ./index/genome.zip -d index/

cp ./index/ncbi_dataset/data/GCF_000195855.1/GCF_000195855.1_ASM19585v1_genomic.fna ./index


# samtool index 
samtools faidx ./index/GCF_000195855.1_ASM19585v1_genomic.fna

#bwa index 
bwa index ./index/GCF_000195855.1_ASM19585v1_genomic.fna -p index/leprae
while IFS=, read sra_id type
do
    #start with data download from NCBI sra
    echo "Downloading : "$sra_id
    prefetch $sra_id -O fastq/rawdata/
    # split the file in and convert fastq
    fasterq-dump -e 4 ./fastq/rawdata/$sra_id -O fastq/rawdata/
    # compress fastq file
    pigz ./fastq/rawdata/*.fastq

    #fastqc for the fastq.gz file
    fastqc ./fastq/rawdata/$sra_id*.fastq.gz -o reports/fastqc_raw

    if [[ "$type" == "SINGLE" ]]
    then
	 echo "SINGLE : "${sra_id}
     fastp -i ./fastq/rawdata/${sra_id}.fastq.gz -o ./fastq/trimmed/${sra_id}_trimmed.fastq.gz -e 20 -l 50 -j ./reports/fastqc_trimmed/${sra_id}_fastp.json -h ./reports/fastqc_trimmed/${sra_id}_fastp.html -w 4

      bwa mem  ./index/leprae ./fastq/trimmed/${sra_id}_trimmed.fastq.gz -t 4 > ./results/aligned_bam/${sra_id}.sam

    elif [[ "$type" == "PAIRED" ]]
    then
    	echo "PAIRED : "${sra_id}
     fastp -i ./fastq/rawdata/${sra_id}_1.fastq.gz -I ./fastq/rawdata/${sra_id}_2.fastq.gz -o ./fastq/trimmed/${sra_id}_1_trimmed.fastq.gz -O ./fastq/trimmed/${sra_id}_2_trimmed.fastq.gz -e 20 -l 50 -h ./reports/fastqc_trimmed/${sra_id}_fastp.html -j ./reports/fastqc_trimmed/${sra_id}_fastp.json -w 4

     bwa mem ./index/leprae ./fastq/trimmed/${sra_id}_1_trimmed.fastq.gz ./fastq/trimmed/${sra_id}_2_trimmed.fastq.gz -t 4 > ./results/aligned_bam/${sra_id}.sam
    fi

    #convert sam to bam
    samtools view -@ 2 -Sb -o ./results/aligned_bam/${sra_id}.bam ./results/aligned_bam/${sra_id}.sam

    # stats count for bam file
    samtools flagstat -@ 2 ./results/aligned_bam/${sra_id}.bam > ./results/aligned_bam/${sra_id}.tsv
    #starting the samtools sort, picard and variant calling
    # sorting
    echo "Sorting : "$sra_id
    samtools sort -@ 2 ./results/aligned_bam/${sra_id}.bam -o ./results/sorted_bam/${sra_id}_sorted.bam -O bam
    # flagstat sorting file
    samtools flagstat -@ 2 ./results/sorted_bam/${sra_id}_sorted.bam > ./reports/samtools_sort/${sra_id}.tsv

    # remove de_duplicates by picard
    picard MarkDuplicates I=./results/sorted_bam/${sra_id}_sorted.bam O=./results/dedup_bam/${sra_id}_dupli.bam M=./reports/samtools_dedup/${sra_id}_marked.txt REMOVE_DUPLICATES=true
    
    # flagstat of picard duplicates
    samtools flagstat -@ 2 ./results/dedup_bam/${sra_id}_dupli.bam > ./reports/dedup_stat/${sra_id}.tsv
    # write a txt for qualimap
    # qualimap multi-bamqc
    
    #mpileup by bcftools

  bcftools mpileup -a "FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR" -O b --threads 2 -o ./results/pileup/${sra_id}_raw.bcf -f ./index/GCF_000195855.1_ASM19585v1_genomic.fna ./results/dedup_bam/${sra_id}_dupli.bam -d 10000

   # varient calling
  bcftools call --ploidy 1 -m -v -o ./results/varient_vcf/${sra_id}_varient.vcf ./results/pileup/${sra_id}_raw.bcf

  #filter the vcf file
  bcftools filter -i 'MQ>30 && INFO/DP>5 && INFO/AD>3' ./results/varient_vcf/${sra_id}_varient.vcf -o ./results/filtered_vcf/${sra_id}_filter.vcf

   #vcf to.gz file
   bgzip ./results/filtered_vcf/${sra_id}_filter.vcf

   # indexing for .gz file
   bcftools index -f -t ./results/filtered_vcf/${sra_id}_filter.vcf.gz



  

done < <(tail -n +2 $1)

    #create a txt file for merge vcf
    ls ./results/filtered_vcf/*_filter.vcf.gz > ./results/anotation/merge.txt
	
      # merge these vcf file by the txt
    bcftools merge -o ./results/anotation/merge.vcf -O v -l ./results/anotation/merge.txt --missing-to-ref

      #normalize genotype by bcftool
    bcftools norm -a -f GCF_000195855.1_ASM19585v1_genomic.fna ./results/anotation/merge.vcf -o ./results/anotation/norm.vcf

      # annotation by snpeff
    snpEff -canon -no-downstream  -no-upstream NC_002677.1 ./results/anotation/norm.vcf -csvStats ./results/anotation/norm_annotate.tsv > ./results/anotation/norm_annotate.vcf

      #extract genotype field
    (bcftools query -l ./results/anotation/norm_annotate.vcf| tr "\n" "\t" &&  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[\t%GT]\n' ./results/anotation/norm_annotate.vcf)> ./results/anotation/genotype_metrixs.tsv

      # extract field by snpshift
     SnpSift extractFields -s "," -e "." ./results/anotation/norm_annotate.vcf CHROM POS REF ALT FILTER ANN[*].GENE ANN[*].GENEID ANN[*].BIOTYPE ANN[*].EFFECT EFF[*].AA EFF[*].EFFECT GEN[*].GT > ./results/anotation/extract_snpshift.tsv

#qualimap text file
for f in $(ls ./results/dedup_bam/*_dupli.bam); do echo -e $(basename ${f})"\t"${f}; done > ./reports/qualimap/qualimap.txt
# multi_bamqc
qualimap multi-bamqc -d ./reports/qualimap/qualimap.txt -outfile ./reports/qualimap/quali_report.pdf -outdir ./reports/qualimap/ -outformat html -r

#merge reports using multiqc
echo "Generating multiqc report"
multiqc ./reports/fastqc_raw/ -o ./reports/fastqc_raw/

#multiqc for trimmed files
multiqc ./reports/fastqc_trimmed/ -o ./reports/fastqc_trimmed/

#multiqc for duplicates remove
multiqc ./reports/dedup_stat/ -o ./reports/dedup_stat/

