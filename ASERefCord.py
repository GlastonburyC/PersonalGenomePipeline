''' This script takes the maternal and paternal coordinate translated BAM files and ASEReadCounterFiles and generates 
    the ASEReadCounter output for each haplotype with the alleles translated onto reference coordinates '''
   
import vcf
import sys
import os

# MATERNAL ASE (GATK output)
ase_mat=open(sys.argv[1],'r')
tmp=[]
snp_mat={}
next(ase_mat)
for line in ase_mat:
	tmp=line.split('\t')
	snp_mat[tmp[2].split(';')[0]]=tmp[1]

# PATERNAL ASE (GATK output)
ase_pat=open(sys.argv[2],'r')
tmp=[]
snp_pat={}
next(ase_pat)
for line in ase_pat:
	tmp=line.split('\t')
	snp_pat[tmp[2].split(';')[0]]=tmp[1]


def merge_two_dicts(x, y):
    '''merge dicts into a new dict as a shallow copy.'''
    z = x.copy()
    z.update(y)
    return z

all_snps=merge_two_dicts(snp_mat,snp_pat)

len(all_snps)
# VCF containing all het SNPs used before ASEReadCounter breaking up into haplotype specific BAMs.
vcf_ref = vcf.Reader(open(sys.argv[3],'r'))

reference_alleles={}
for record in vcf_ref:
	if record.ID.split(';')[0] in all_snps:
		reference_alleles[record.ID.split(';')[0]] = record.REF,str(record.ALT[0]),record.POS

len(reference_alleles)
'''
	ASEReadCounter outputs indels as A/A when in fact it should be AAGG / A for example. 
        Correct alleles with full nucleotide string. (GATK not fixed as of 08/04/16)
'''
ase_mat=open(sys.argv[1],'r')

vcf_mat = vcf.Reader(open(sys.argv[4],'r'))
vcf_pat = vcf.Reader(open(sys.argv[5],'r'))

ase_mat_out=open(sys.argv[6],'w')

maternal_store={}
for record in vcf_mat:
	A1 = record.genotype(sys.argv[10])['GT'].split('|')[0]
	A2 = record.genotype(sys.argv[10])['GT'].split('|')[1]
	maternal_store[record.ID] = record.REF,record.ALT[0],A1,str(A2[0])

paternal_store={}
for record2 in vcf_pat:
	A1 = record2.genotype(sys.argv[10])['GT'].split('|')[0]
	A2 = record2.genotype(sys.argv[10])['GT'].split('|')[1]
	paternal_store[record2.ID] = record2.REF,record2.ALT[0],A1,str(A2[0])

header=next(ase_mat)
ase_mat_out.write(header)
for line in (line.strip().split() for line in ase_mat):
	if line[2] in maternal_store:
			ase_mat_out.write('\t'.join(line[0:3])+'\t'+str(maternal_store[line[2]][0])+"\t"+str(maternal_store[line[2]][1])+'\t'+'\t'.join(line[5:])+'\n')

ase_mat_out.close()

ase_pat=open(sys.argv[2],'r')
ase_pat_out=open(sys.argv[7],'w')
header=next(ase_pat)
ase_pat_out.write(header)

for line in (line.strip().split() for line in ase_pat):
	if line[2] in paternal_store:
			ase_pat_out.write('\t'.join(line[:3])+'\t'+str(paternal_store[line[2]][0])+"\t"+str(paternal_store[line[2]][1])+'\t'+'\t'.join(line[5:])+'\n')

ase_pat_out.close()

os.system('mv %s %s' % (sys.argv[7], sys.argv[2]))
os.system('mv %s %s' % (sys.argv[6], sys.argv[1]))

######### convert to reference alleles and reference position. ###########

ase_mat=open(sys.argv[1],'r')

out_ase=open(sys.argv[8],'w')

header=next(ase_mat)
out_ase.write(header)
for mat in ase_mat:
	if mat.split('\t')[2].split(';')[0] in reference_alleles:
		if mat.split('\t')[3] != reference_alleles[mat.split('\t')[2].split(';')[0]][0]:
			out_ase.write(mat.split('\t')[0]+'\t'+str(reference_alleles[mat.split('\t')[2].split(';')[0]][2])+'\t'+mat.split('\t')[2].split(';')[0]+'\t'+str(reference_alleles[mat.split('\t')[2].split(';')[0]][0])+'\t'+str(reference_alleles[mat.split('\t')[2].split(';')[0]][1])+'\t'+mat.split('\t')[6]+'\t'+mat.split('\t')[5]+'\t'+'\t'.join(mat.split('\t')[7:]))
		else:
			out_ase.write(mat.split('\t')[0]+'\t'+str(reference_alleles[mat.split('\t')[2].split(';')[0]][2])+'\t'+mat.split('\t')[2].split(';')[0]+'\t'+str(reference_alleles[mat.split('\t')[2].split(';')[0]][0])+'\t'+str(reference_alleles[mat.split('\t')[2].split(';')[0]][1])+'\t'+mat.split('\t')[5]+'\t'+mat.split('\t')[6]+'\t'+'\t'.join(mat.split('\t')[7:]))

out_ase.close()

ase_pat=open(sys.argv[2],'r')

out_ase=open(sys.argv[9],'w')

header=next(ase_pat)
out_ase.write(header)
for pat in ase_pat:
	if pat.split('\t')[2].split(';')[0] in reference_alleles:
		if pat.split('\t')[3] != reference_alleles[pat.split('\t')[2].split(';')[0]][0]:
			out_ase.write(pat.split('\t')[0]+'\t'+str(reference_alleles[pat.split('\t')[2].split(';')[0]][2])+'\t'+pat.split('\t')[2].split(';')[0]+'\t'+str(reference_alleles[pat.split('\t')[2].split(';')[0]][0])+'\t'+str(reference_alleles[pat.split('\t')[2].split(';')[0]][1])+'\t'+pat.split('\t')[6]+'\t'+pat.split('\t')[5]+'\t'+'\t'.join(pat.split('\t')[7:]))
		else:
			out_ase.write(pat.split('\t')[0]+'\t'+str(reference_alleles[pat.split('\t')[2].split(';')[0]][2])+'\t'+pat.split('\t')[2].split(';')[0]+'\t'+str(reference_alleles[pat.split('\t')[2].split(';')[0]][0])+'\t'+str(reference_alleles[pat.split('\t')[2].split(';')[0]][1])+'\t'+pat.split('\t')[5]+'\t'+pat.split('\t')[6]+'\t'+'\t'.join(pat.split('\t')[7:]))

out_ase.close()
