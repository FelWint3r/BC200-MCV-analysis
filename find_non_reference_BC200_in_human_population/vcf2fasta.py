of=open('shapeit5-phased-callset_final-vcf.phased.vcf','r').readlines()
for i in of:
 if not (i.startswith('#')):
  temp=i.split('\t')
  chrno=temp[0]
  position=temp[1]
  seq=temp[4]
  title='>'+chrno+'_'+position+'_'
  if(len(seq)>80 and len(seq)<300):
   print(title)
   print(seq)