#halBranchMutations /mnt/sunLab/gaopu/mcv/bc200history/test_cactus/run_halbranchmutations/8primates.hal  Homo_sapiens --refFile Homo_sapiens_ins.bed --parentFile Homo_sapiens_del.bed
python3 bed2seq.py > insseq.fasta &&
RepeatMasker -lib /mnt/sunLab/gaopu/mcv/bc200history/test_cactus/consensus_alu.fasta -engine crossmatch -parallel 20 insseq.fasta &&
grep BC200 insseq.fasta.out | awk '{print $5,$6,$7,$9}' > repeatmasker_output_BC200.txt &&
python3 ins_site2sequence.py &&
cat BC200_gene_consensus.fasta ins_sequence.fasta > ins_sequence_with_BC200gene.fasta &&
clustalw2 -INFILE=ins_sequence_with_BC200gene.fasta -OUTORDER=INPUT -OUTPUT=FASTA -OUTFILE=ins_sequence_with_BC200gene_aligned.fasta &&
python3 genome_coordinate_correction.py &&
python3 diagnostic_sites_count_with_tsd_finder.py > ins_with_at_least10_diagnostic_site.txt &&
cat *water* >total_tsd_file.txt
