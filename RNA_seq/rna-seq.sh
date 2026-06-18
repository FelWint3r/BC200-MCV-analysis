# Quality control

java -jar /mnt/sunLab/TJX/software/Trimmomatic/trimmomatic-0.40.jar 
SE -phred33 
SRR3142259.fastq 
SRR3142259.clean.fastq.gz 
ILLUMINACLIP:/mnt/sunLab/TJX/software/Trimmomatic/adapters/TruSeq3-SE.fa:2:30:10 
LEADING:3 
TRAILING:3 
SLIDINGWINDOW:4:15 
MINLEN:36

# Alignment

bowtie2 
-x ../bowtie2_db/db 
-U SRR3142259.clean.fastq.gz 
-S SRR3142259.sam

# SAM to BAM conversion

samtools view -bS SRR3142259.sam \

> SRR3142259.bam

# Retain alignments with MAPQ > 10

samtools view -h -q 10 SRR3142259.bam \

> SRR3142259.mapq10.sam

# Convert filtered SAM to BAM

samtools view -bS SRR3142259.mapq10.sam \

> SRR3142259.mapq10.bam

# Sort BAM

samtools sort 
SRR3142259.mapq10.bam 
-o SRR3142259.mapq10.sorted.bam

# Index BAM

samtools index 
SRR3142259.mapq10.sorted.bam

# Gene expression quantification

featureCounts 
-T 10 
-M 
-O 
--fraction 
-t exon 
-g gene 
-a /mnt/sunLab/tanshuxin/BC200/human/all_human/t2t.gff 
-o SRR3142259.featureCounts.txt 
SRR3142259.mapq10.sorted.bam










