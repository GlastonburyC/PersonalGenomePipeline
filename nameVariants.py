import vcf
import sys
vcf_ref = vcf.Reader(open(sys.argv[1]+'/'+sys.argv[2],'r'))
vcf_mat_out = vcf.Writer(open(sys.argv[1]+'/'+sys.argv[3],'w'),vcf_ref)

reference_alleles={}
for record in vcf_ref:
    if record.ID == None:
    	record.ID = str(record.CHROM)+':'+str(record.POS)
    	vcf_mat_out.write_record(record)
    else:
    	vcf_mat_out.write_record(record)
