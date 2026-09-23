from pyfaidx import Fasta

genome="/mnt/sunLab/gaopu/mcv/bc200history/test_cactus/homo_sapiens/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"

fasta=Fasta(genome)

of=open('Homo_sapiens_ins.bed','r').readlines()
for i in of:
 if not (i.startswith('#')):
  temp=i.split('\t')
  chrno=temp[0]
  l_site=int(temp[1])
  r_site=int(temp[2])
  svtype=temp[3]
  length=r_site-l_site+1
  if(svtype=='I' and length>=135 and length<=250):
   seq=fasta[chrno][l_site-1:r_site]
   print('>'+chrno+'*'+str(l_site)+'*'+str(r_site))
   print(seq)




