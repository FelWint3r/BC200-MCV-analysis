# -*- coding: utf-8 -*-
import warnings
warnings.filterwarnings("ignore")
from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline
from pyfaidx import Fasta
import subprocess
import os

fasta = Fasta("/mnt/sunLab/gaopu/mcv/bc200history/test_cactus/samealu/8_live_pri/homo_sapiens/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna")


def get_nth_base_position_in_aligned(aligned_seq, n):
    count = 0
    for i, base in enumerate(aligned_seq):
        if base.upper() in "ATGC":
            count += 1
            if count == n:
                return i
    return -1

BC200="GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCTCTTAGGGAGGCAGAGGCGGGAGGATAGCTTGAGCCCAGGAGTTCGAGACCTGCCTGGGCAATATAGCGAGACCCCGTTCTCCACAAAAAGGAAAAAAAAAAAAAAAA"

of=open('ins_info.txt','r').readlines()

for i in of:
 temp=i.split('*')
 strand=temp[3]
 l_site=int(temp[1])
 r_site=int(temp[2])
 chrno=temp[0][1:]
 title=temp[0][1:]+'*'+temp[1]+'*'+temp[2]+'*'+temp[3]
 seq=temp[4].upper()
 old_seq=temp[4]
 temp_file_name='temp'+'*'+title+'.fasta'
 final_file_name='final'+'*'+title+'.fasta'
 tf=open(temp_file_name,'w')
 tf.write('>BC200'+'\n')
 tf.write(BC200+'\n')
 tf.write('>'+title+'\n')
 tf.write(seq+'\n')
 tf.close()
#run clustalw
 clustalw_cline = ClustalwCommandline("clustalw2", infile=temp_file_name, outfile=final_file_name, output="fasta")
 clustalw_cline()
#process clustalw result
 for item in SeqIO.parse(final_file_name, "fasta"):
  if(item.id=='BC200'):
   seq=str(item.seq)
   p1=get_nth_base_position_in_aligned(item.seq,36)
   p2=get_nth_base_position_in_aligned(item.seq,40)
   p3=get_nth_base_position_in_aligned(item.seq,48)
   p4=get_nth_base_position_in_aligned(item.seq,62)
   p5=get_nth_base_position_in_aligned(item.seq,87)
   p6=get_nth_base_position_in_aligned(item.seq,98)
   p7=get_nth_base_position_in_aligned(item.seq,113)
   p8=get_nth_base_position_in_aligned(item.seq,118)
   p9=get_nth_base_position_in_aligned(item.seq,120)
   p10=get_nth_base_position_in_aligned(item.seq,126)
   p11=get_nth_base_position_in_aligned(item.seq,127)
  elif(item.id!='BC200'):
   seq=str(item.seq)
   score=0
   ds=item.seq[p1]+item.seq[p2]+item.seq[p3]+item.seq[p4]+item.seq[p5]+item.seq[p6]+item.seq[p7]+item.seq[p8]+item.seq[p9]+item.seq[p10]+item.seq[p11]
   if(ds[0]=='T'):
    score=score+1
   if(ds[1]=='A'):
    score=score+1
   if(ds[2]=='A'):
    score=score+1
   if(ds[3]=='A'):
    score=score+1
   if(ds[4]=='T'):
    score=score+1
   if(ds[5]=='T'):
    score=score+1
   if(ds[6]=='T'):
    score=score+1
   if(ds[7]=='C'):
    score=score+1
   if(ds[8]=='G'):
    score=score+1
   if(ds[9]=='G'):
    score=score+1
   if(ds[10]=='G'):
    score=score+1
    
    
   if(score>=10):
    if(strand=='C'):
     up_seq=str(fasta[chrno][r_site:r_site+30]).upper()
     down_seq=str(fasta[chrno][l_site-31:l_site+9]).upper()
    elif(strand=='+'):
     up_seq=str(fasta[chrno][l_site-31:l_site-1]).upper()
     down_seq=str(fasta[chrno][r_site-9:r_site+31]).upper()
     
    up_name=title+'*up.fasta'
    down_name=title+'*down.fasta'
    out_name=title+'*water.txt'
    
    open_up=open(up_name,'w')
    open_down=open(down_name,'w')   
    
    open_up.write('>'+title+'*up'+'\n'+up_seq)
    open_down.write('>'+title+'*down'+'\n'+down_seq)
    
    open_up.close()
    open_down.close()
    
    cmd='water'+' '+'-asequence'+' '+up_name+' '+'-bsequence'+' '+down_name+' '+'-stdout'+' '+'-auto'+' '+'-outfile'+' '+out_name

    
    os.system(cmd)
    open_water_result=open(out_name,'r').readlines()
    match=''
    for record in open_water_result:     
     if('|' in record):
      match=match+(record[21:-1])

    if(len(match)>=10):    
     for y in range(len(match) - 9):
      if match[y:y+10].count('|') >= 7:
       l_site_final=temp[1]
       r_site_final=temp[2]
       print(temp[0][1:]+'\t'+l_site_final+'\t'+r_site_final+'\t'+temp[3]+'\t'+' '+'\t'+str(score)+'\t'+old_seq[:-1])

       break
    elif(len(match)<10):
     if(match.count('|'))>=7:
      l_site_final=temp[1]
      r_site_final=temp[2]
      print(temp[0][1:]+'\t'+l_site_final+'\t'+r_site_final+'\t'+temp[3]+'\t'+' '+'\t'+str(score)+'\t'+old_seq[:-1])

        
      
    
#os.system('rm -fr *up* *down*')
#os.system('rm -fr *final* *temp*')
#os.system('rm -fr *water*')
   
