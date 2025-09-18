#!/usr/bin/bash

# this script is to test multiple levels of 
# dada2-min-abundance

# set the working directory
cd ~/projects/abbott-carbohydrate-obesity-project2/microbiome_analysis

# if the output file asvs.fasta does not exist then run this part of the analysis.

if [ ! -f dada2_asv_1/asvs.fasta ]; then
    ./asv-finder derep2asv --method dada2 --input-dir fastq_filtered/ --output-dir dada2_asv_1 --threads 20 --dada2-min-abundance 1 --dada2-adaptive-screening
fi

if [ ! -f dada2_asv_2/asvs.fasta ]; then
    ./asv-finder derep2asv --method dada2 --input-dir fastq_filtered/ --output-dir dada2_asv_2 --threads 20 --dada2-min-abundance 2 --dada2-adaptive-screening
fi

if [ ! -f dada2_asv_3/asvs.fasta ]; then
    ./asv-finder derep2asv --method dada2 --input-dir fastq_filtered/ --output-dir dada2_asv_3 --threads 20 --dada2-min-abundance 3 --dada2-adaptive-screening
fi

if [ ! -f dada2_asv_4/asvs.fasta ]; then
    ./asv-finder derep2asv --method dada2 --input-dir fastq_filtered/ --output-dir dada2_asv_4 --threads 20 --dada2-min-abundance 4 --dada2-adaptive-screening
fi

if [ ! -f dada2_asv_5/asvs.fasta ]; then
    ./asv-finder derep2asv --method dada2 --input-dir fastq_filtered/ --output-dir dada2_asv_5 --threads 20 --dada2-min-abundance 5 --dada2-adaptive-screening
fi

if [ ! -f dada2_asv_6/asvs.fasta ]; then
    ./asv-finder derep2asv --method dada2 --input-dir fastq_filtered/ --output-dir dada2_asv_6 --threads 20 --dada2-min-abundance 6 --dada2-adaptive-screening
fi

if [ ! -f dada2_asv_7/asvs.fasta ]; then
    ./asv-finder derep2asv --method dada2 --input-dir fastq_filtered/ --output-dir dada2_asv_7 --threads 20 --dada2-min-abundance 7 --dada2-adaptive-screening
fi
 
# now if the output file strobelca_taxonomy.tsv does not exist then run this part of the analysis.

if [ ! -f dada2_asv_1/strobelca_taxonomy.tsv ]; then
    ./asv-finder taxonomy --method strobelca --input-fasta dada2_asv_1/asvs.fasta --output dada2_asv_1/strobelca_taxonomy --threads 20 --database references/SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz
fi

if [ ! -f dada2_asv_2/strobelca_taxonomy.tsv ]; then
    ./asv-finder taxonomy --method strobelca --input-fasta dada2_asv_2/asvs.fasta --output dada2_asv_2/strobelca_taxonomy --threads 20 --database references/SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz
fi

if [ ! -f dada2_asv_3/strobelca_taxonomy.tsv ]; then
    ./asv-finder taxonomy --method strobelca --input-fasta dada2_asv_3/asvs.fasta --output dada2_asv_3/strobelca_taxonomy --threads 20 --database references/SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz
fi

if [ ! -f dada2_asv_4/strobelca_taxonomy.tsv ]; then
    ./asv-finder taxonomy --method strobelca --input-fasta dada2_asv_4/asvs.fasta --output dada2_asv_4/strobelca_taxonomy --threads 20 --database references/SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz
fi

if [ ! -f dada2_asv_5/strobelca_taxonomy.tsv ]; then
    ./asv-finder taxonomy --method strobelca --input-fasta dada2_asv_5/asvs.fasta --output dada2_asv_5/strobelca_taxonomy --threads 20 --database references/SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz
fi

if [ ! -f dada2_asv_6/strobelca_taxonomy.tsv ]; then
    ./asv-finder taxonomy --method strobelca --input-fasta dada2_asv_6/asvs.fasta --output dada2_asv_6/strobelca_taxonomy --threads 20 --database references/SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz
fi

if [ ! -f dada2_asv_7/strobelca_taxonomy.tsv ]; then
    ./asv-finder taxonomy --method strobelca --input-fasta dada2_asv_7/asvs.fasta --output dada2_asv_7/strobelca_taxonomy --threads 20 --database references/SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz
fi
