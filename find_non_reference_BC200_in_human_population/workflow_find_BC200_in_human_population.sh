#gunzip shapeit5-phased-callset_final-vcf.phased.vcf.gz
python3 vcf2fasta.py >vcf.fasta &&
makeblastdb  -dbtype nucl -out db -in vcf.fasta &&
blastn -task blastn -query Human_BC200_gene.fasta  -out find_non_reference_BC200.txt -db db -evalue 1e-30