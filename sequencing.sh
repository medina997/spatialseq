#! /bin/bash
SECONDS= 0

#path
SAMPLES="/media/hillierlab/6b8db916-205e-4bc2-bc7a-ac183d8928f1/data/Test/macaque4_samples"
TRIM="/home/hillierlab/Downloads/Trimmomatic-0.39/trimmomatic-0.39.jar"
REF="/media/hillierlab/6b8db916-205e-4bc2-bc7a-ac183d8928f1/data/Test/Macaque_reference_genome/Macaca_fascicularis_6.0.dna.toplevel.fa"
GTF="/media/hillierlab/6b8db916-205e-4bc2-bc7a-ac183d8928f1/data/Test/Macaque_reference_genome/Macaca_fascicularis_6.0.107.gtf"

cd $SAMPLES

#Quality analysis
#fastqc 2-A2_S14_R1_001.fastq.gz

#trimming the reads
java -jar $TRIM SE -threads 10 -phred33 2-A3_S15_R1_001.fastq.gz clean3.fastq CROP:75 MINLEN:60

#Quality analysis after trimming
fastqc clean3.fastq 

#STAR alignment
#Creating the genome index
#only need to do once
#create folder where the indices are stored
#mkdir STAR_alignment_sample2
#STAR --runMode genomeGenerate --runThreadN 2 --genomeDir #STAR_alignment_sample2 --genomeFastaFiles $REF --sjdbGTFfile $GFF --#sjdbOverhang 50 --limitGenomeGenerateRAM 35914453718

#STAR alignment
STAR --quantMode GeneCounts --genomeDir /media/hillierlab/6b8db916-205e-4bc2-bc7a-ac183d8928f1/data/Test/Macaque_reference_genome/STAR_macdna_6.0.107.gtf --runThreadN 10 --readFilesIn /media/hillierlab/6b8db916-205e-4bc2-bc7a-ac183d8928f1/data/Test/macaque4_samples/clean3.fastq --outFileNamePrefix wt3_ --outSAMtype BAM SortedByCoordinate

#with multiqc check the alignment
#multiqc .

#FeatureCounts
mkdir FeatureCount_sample3_gtf
featureCounts -T 1 -g gene_id --fracOverlap 0.25 -a $GTF -o FeatureCount_sample3_gtf/readcounts.txt *.bam
