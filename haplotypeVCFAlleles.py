import vcf 
import os
import sys

vcf_mat = vcf.Reader(open(sys.argv[1]+'/maternal/'+sys.argv[2],'r'))

vcf_mat_out = vcf.Writer(open(sys.argv[1]+'/maternal/'+sys.argv[2]+'.maternal2.vcf','w'),vcf_mat)

for record in vcf_mat:
	if record.genotype(sys.argv[1])['GT'].split('|')[1] == '1':
		tmp=record
		A1=tmp.REF
		A2=tmp.ALT
		tmp.REF=str(A2[0])
		tmp.ALT=[vcf.model._Substitution(A1)]
		vcf_mat_out.write_record(tmp)
	else:
		vcf_mat_out.write_record(record)

vcf_mat_out.close()
#os.system('%s 131.maternal.vcf' % sys.argv[2])

import vcf 

vcf_pat = vcf.Reader(open(sys.argv[1]+'/paternal/'+sys.argv[3],'r'))
vcf_pat_out = vcf.Writer(open(sys.argv[1]+'/paternal/'+sys.argv[3]+'.paternal2.vcf','w'),vcf_pat)

for record in vcf_pat:
	if record.genotype(sys.argv[1])['GT'].split('|')[0] == '1':
		tmp=record
		A1=tmp.REF
		A2=tmp.ALT
		tmp.REF=str(A2[0])
		tmp.ALT=[vcf.model._Substitution(A1)]
		vcf_pat_out.write_record(tmp)
	else:
		vcf_pat_out.write_record(record)

vcf_pat_out.close()
#os.system('%s 131.paternal.vcf' % sys.argv[2])
