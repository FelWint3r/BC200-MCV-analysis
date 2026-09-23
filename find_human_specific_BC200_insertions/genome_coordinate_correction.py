from Bio import SeqIO

from pyfaidx import Fasta
fasta = Fasta("/mnt/sunLab/gaopu/mcv/bc200history/test_cactus/samealu/8_live_pri/homo_sapiens/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna")
def reverse_complement(dna_sequence):
  
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C','a':'t','t':'a','g':'c','c':'g'}
    

    complement_seq = [complement[base] for base in dna_sequence]
    
    return ''.join(complement_seq)[::-1]
count=0
f_reviewed=open('ins_info.txt','w')

for record in SeqIO.parse("ins_sequence_with_BC200gene_aligned.fasta", "fasta"):
 seq = str(record.seq)
 title=record.id
 if(title=='BC200#BC200_gene_consensus_sequence'):
  for i in seq:
   if(i=='-'):
    count=count+1
   else:
    break
 elif(title!='BC200#BC200_gene_consensus_sequence'):
  temp=title.split('*')
  chrno=temp[0]
  l_site=int(temp[1])
  r_site=int(temp[2])
  strand=temp[3]
  duoyu=count-seq[:count].count('-')
  
  if(strand=='C'):
  
   new_seq=reverse_complement(str(fasta[chrno][l_site-1:r_site-duoyu]))
   
   new_title='>'+chrno+'*'+str(l_site)+'*'+str(r_site-duoyu)+'*'+'C'
   
   
   
   f_reviewed.write(new_title+'*'+new_seq+'\n')
    
  elif(strand=='+'):
  
   new_seq=str(fasta[chrno][l_site+duoyu-1:r_site])
   
   new_title='>'+chrno+'*'+str(l_site+duoyu)+'*'+str(r_site)+'*'+'+'
   
   
   
   f_reviewed.write(new_title+'*'+new_seq+'\n')
    
f_reviewed.close()