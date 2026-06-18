import re
from pyfaidx import Fasta
fasta = Fasta("/mnt/sunLab/gaopu/mcv/bc200history/test_cactus/samealu/8_live_pri/homo_sapiens/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna")
def reverse_complement(dna_sequence):
  
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C','a':'t','t':'a','g':'c','c':'g'}
    

    complement_seq = [complement[base] for base in dna_sequence]
    
    return ''.join(complement_seq)[::-1]
of=open('repeatmasker_output_BC200.txt','r').readlines()
f_w1=open('ins_sequence.fasta','w')
for i in of:
 temp=i.split(' ')
 chrno=temp[0].split('*')[0]
 start=int(temp[0].split('*')[1])
 end=int(temp[0].split('*')[2])
 seq_s=int(temp[1])
 seq_e=int(temp[2])
 strand=temp[3][:-1]
 l_site=start+seq_s-1
 r_site=start+seq_e-1
 if(strand=='C'):
   seq=reverse_complement(str(fasta[chrno][l_site:r_site]))
   f_w1.write('>'+chrno+'*'+str(l_site)+'*'+str(r_site)+'*'+strand+'*'+'\n')
   f_w1.write(seq+'\n') 
   print(l_site,r_site)
 if(strand=='+'):
   seq=str(fasta[chrno][l_site:r_site])
   f_w1.write('>'+chrno+'*'+str(l_site)+'*'+str(r_site)+'*'+strand+'*'+'\n')
   f_w1.write(seq+'\n')
f_w1.close()
