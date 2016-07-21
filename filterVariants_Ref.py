import vcf
import sys
import os

ID=str(sys.argv[1])
vcf_reader = vcf.Reader(open(sys.argv[1]+'/'+sys.argv[1]+'.vcf.gz','r'))
vcf_writer = vcf.Writer(open(sys.argv[1]+'/'+sys.argv[1]+'.hets.phased.vcf','w'),vcf_reader)
for record in vcf_reader:
  if len(record.get_hets()) != 0:
  	if record.genotype(ID)['GT'][1]=='|':
  		if sorted(record.genotype(ID)['PL'])[1] >= 30:
  			vcf_writer.write_record(record)
vcf_writer.close()
