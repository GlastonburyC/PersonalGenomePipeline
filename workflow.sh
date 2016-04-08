# Maternal
STAR/bin/Linux_x86_64/STAR --runThreadN 8 --runMode genomeGenerate --genomeDir maternal/ --genomeFastaFiles *_131_maternal.fa

#Paternal 
STAR/bin/Linux_x86_64/STAR --runThreadN 8 --runMode genomeGenerate --genomeDir paternal/ --genomeFastaFiles *_131_paternal.fa

#Align paternal
STAR/bin/Linux_x86_64/STAR --runThreadN 8 --runMode alignReads --readFilesIn EB131/F226_new_120618_7_sorted.bam.f1_val_1.fq EB131/F226_new_120618_7_sorted.bam.f2_val_2.fq --genomeDir paternal --outFilterMultimapNmax 30 --sjdbOverhang 48 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattrIHstart 0 --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --chimSegmentMin 15 --outMultimapperOrder Random --outFilterMismatchNmax 4 --outSAMunmapped Within --outSAMattributes NH HI NM MD

#Align maternal
STAR/bin/Linux_x86_64/STAR --runThreadN 8 --runMode alignReads --readFilesIn EB131/F226_new_120618_7_sorted.bam.f1_val_1.fq EB131/F226_new_120618_7_sorted.bam.f2_val_2.fq --genomeDir maternal --outFilterMultimapNmax 30 --sjdbOverhang 48 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattrIHstart 0 --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --chimSegmentMin 15 --outMultimapperOrder Random --outFilterMismatchNmax 4 --outSAMunmapped Within --outSAMattributes NH HI NM MD

########################################

vcf-sort -c 131.hets.GATK.vcf > 131.hets.GATK.sorted.vcf

#########################################

#rm paternal.edit.chain
#rm maternal.edit.chain
#rm all_paternal.fa
#rm all_maternal.fa
	#statements

# run BASH script to concatenate reference genomes, and rename chain 1_maternal > 1

cp paternal.chain paternal.edit.chain
cp maternal.chain maternal.edit.chain
for i in {1..22}
do
	sed -i 's/'"$i"'_maternal/'"$i"'/g' maternal.edit.chain
	sed -i 's/'"$i"'_paternal/'"$i"'/g' paternal.edit.chain
	cat "$i"_131_paternal.fa >> all_paternal.fa
	cat "$i"_131_maternal.fa >> all_maternal.fa

done

sed -i 's/X_paternal/X/g' paternal.edit.chain
sed -i 's/Y_paternal/Y/g' paternal.edit.chain
sed -i 's/MT_paternal/MT/g' paternal.edit.chain

sed -i 's/X_maternal/X/g' maternal.edit.chain
sed -i 's/Y_maternal/Y/g' maternal.edit.chain
sed -i 's/MT_maternal/MT/g' maternal.edit.chain

cat X_131_paternal.fa >> all_paternal.fa
cat Y_131_paternal.fa >> all_paternal.fa
cat MT_131_paternal.fa >> all_paternal.fa

cat X_131_maternal.fa >> all_maternal.fa
cat Y_131_maternal.fa >> all_maternal.fa
cat MT_131_maternal.fa >> all_maternal.fa

#################################


# Remap VCF to parental genome coordinates using modified CrossMap (fixedBugs)
sudo python mod.py vcf maternal.edit.chain 131.hets.GATK.sorted.renamed.vcf maternal.renamed.fa 131.maternal.vcf

sudo python mod.py vcf paternal.edit.chain 131.hets.GATK.sorted.renamed.vcf paternal.renamed.fa 131.paternal.vcf


import vcf 
import os
vcf_mat = vcf.Reader(open('131.maternal.vcf','r'))
vcf_mat_out = vcf.Writer(open('131.maternal2.vcf','w'),vcf_mat)

for record in vcf_mat:
	if record.genotype('131')['GT'].split('|')[1] == '1':
		tmp=record
		A1=tmp.REF
		A2=tmp.ALT
		tmp.REF=str(A2[0])
		tmp.ALT=[vcf.model._Substitution(A1)]
		vcf_mat_out.write_record(tmp)
	else:
		vcf_mat_out.write_record(record)

vcf_mat_out.close()
os.system('mv 131.maternal2.vcf 131.maternal.vcf')

import vcf 

vcf_pat = vcf.Reader(open('131.paternal.vcf','r'))
vcf_pat_out = vcf.Writer(open('131.paternal2.vcf','w'),vcf_pat)

for record in vcf_pat:
	if record.genotype('131')['GT'].split('|')[0] == '1':
		tmp=record
		A1=tmp.REF
		A2=tmp.ALT
		tmp.REF=str(A2[0])
		tmp.ALT=[vcf.model._Substitution(A1)]
		vcf_pat_out.write_record(tmp)
	else:
		vcf_pat_out.write_record(record)

vcf_pat_out.close()
os.system('mv 131.paternal2.vcf 131.paternal.vcf')





/usr/lib/jvm/java-8-oracle/bin/java -jar picard-tools-2.1.1/picard.jar CreateSequenceDictionary R=paternal.renamed.fa O=paternal.renamed.dict
/usr/lib/jvm/java-8-oracle/bin/java -jar picard-tools-2.1.1/picard.jar CreateSequenceDictionary R=maternal.renamed.fa O=maternal.renamed.dict

/usr/lib/jvm/java-8-oracle/bin/java -jar picard-tools-2.1.1/picard.jar AddOrReplaceReadGroups I=EB131/consensus.mat.filtered.sorted.bam O=EB131/consensus.mat.filtered.sorted.readGroup.bam  RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
/usr/lib/jvm/java-8-oracle/bin/java -jar picard-tools-2.1.1/picard.jar AddOrReplaceReadGroups I=EB131/consensus.pat.filtered.sorted.bam O=EB131/consensus.pat.filtered.sorted.readGroup.bam  RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20

samtools view -h EB131/consensus.mat.filtered.sorted.bam | sed -e 's/_maternal//g' >> EB131/consensus.mat.filtered.sorted2.bam
samtools view -h EB131/consensus.pat.filtered.sorted.readGroup2.bam | sed -e 's/_maternal//g' >> EB131/consensus.pat.filtered.sorted2.bam

/usr/lib/jvm/java-8-oracle/bin/java -jar picard-tools-2.1.1/picard.jar ReorderSam I=EB131/consensus.mat.filtered.sorted.readGroup2.bam O=EB131/consensus.mat.filtered.sorted.readGroup2.sorted.bam R=maternal.renamed.fa ALLOW_CONTIG_LENGTH_DISCORDANCE=true
/usr/lib/jvm/java-8-oracle/bin/java -jar picard-tools-2.1.1/picard.jar ReorderSam I=EB131/consensus.pat.filtered.sorted2.bam O=EB131/consensus.pat.filtered.sorted.readGroup2.sorted.bam R=paternal.renamed.fa ALLOW_CONTIG_LENGTH_DISCORDANCE=true

java -jar GenomeAnalysisTK.jar -R maternal.renamed.fa -T ASEReadCounter -o ASE.mat.csv -I EB131/consensus.mat.filtered.sorted.readGroup2.sorted.bam -sites 131.maternal.vcf -dels -U ALLOW_N_CIGAR_READS -S SILENT
java -jar GenomeAnalysisTK.jar -R paternal.renamed.fa -T ASEReadCounter -o ASE.pat.csv -I EB131/consensus.pat.filtered.sorted.readGroup2.sorted.bam -sites 131.paternal.vcf -dels -U ALLOW_N_CIGAR_READS -S SILENT


### Currently, the ASE output is relative to the haplotype reference used (Maternal reference or paternal reference) i.e. REF in ASE.mat.csv = maternal allele etc.
### To assess REF Bias, and for easier interpretation, it's better to have ASE in terms of universal REF / ALT. i.e. all adding up all the REF/ALT counts for each haplotype, consolidated into
### a single ASE file.

import vcf

ase_mat=open('ASE.mat.csv','r')
tmp=[]
snp_mat={}
next(ase_mat)
for line in ase_mat:
	tmp=line.split('\t')
	snp_mat[tmp[2].split(';')[0]]=tmp[1]

ase_pat=open('ASE.pat.csv','r')
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


vcf_ref = vcf.Reader(open('131.hets.GATK.sorted.renamed.vcf','r'))

reference_alleles={}
for record in vcf_ref:
	if record.ID.split(';')[0] in all_snps:
		reference_alleles[record.ID.split(';')[0]] = record.REF,str(record.ALT[0]),record.POS

len(reference_alleles)

######### ASEReadCounter outputs indels as A/A when in fact it should be AAGG / A for example. Correct alleles with full nucleotide string ###########

ase_mat=open('ASE.mat.csv','r')

vcf_mat = vcf.Reader(open('131.maternal.vcf','r'))
vcf_pat = vcf.Reader(open('131.paternal.vcf','r'))

ase_mat_out=open('ASE.mat.alleles.csv','w')

maternal_store={}
for record in vcf_mat:
	A1 = record.genotype('131')['GT'].split('|')[0]
	A2 = record.genotype('131')['GT'].split('|')[1]
	maternal_store[record.ID] = record.REF,record.ALT[0],A1,str(A2[0])

paternal_store={}
for record2 in vcf_pat:
	A1 = record2.genotype('131')['GT'].split('|')[0]
	A2 = record2.genotype('131')['GT'].split('|')[1]
	paternal_store[record2.ID] = record2.REF,record2.ALT[0],A1,str(A2[0])

header=next(ase_mat)
ase_mat_out.write(header)
for line in (line.strip().split() for line in ase_mat):
	if line[2] in maternal_store:
			ase_mat_out.write('\t'.join(line[0:3])+'\t'+str(maternal_store[line[2]][0])+"\t"+str(maternal_store[line[2]][1])+'\t'+'\t'.join(line[5:])+'\n')

ase_mat_out.close()

ase_pat=open('ASE.pat.csv','r')
ase_pat_out=open('ASE.pat.alleles.csv','w')
header=next(ase_pat)
ase_pat_out.write(header)

for line in (line.strip().split() for line in ase_pat):
	if line[2] in paternal_store:
			ase_pat_out.write('\t'.join(line[:3])+'\t'+str(paternal_store[line[2]][0])+"\t"+str(paternal_store[line[2]][1])+'\t'+'\t'.join(line[5:])+'\n')

ase_pat_out.close()

os.system('mv ASE.pat.alleles.csv ASE.pat.csv')
os.system('mv ASE.mat.alleles.csv ASE.mat.csv')

######### convert to reference alleles and reference position. ###########

ase_mat=open('ASE.mat.alleles.csv','r')

out_ase=open('ASE.mat.refAllele.csv','w')

header=next(ase_mat)
out_ase.write(header)
for mat in ase_mat:
	if mat.split('\t')[2].split(';')[0] in reference_alleles:
		if mat.split('\t')[3] != reference_alleles[mat.split('\t')[2].split(';')[0]][0]:
			out_ase.write(mat.split('\t')[0]+'\t'+str(reference_alleles[mat.split('\t')[2].split(';')[0]][2])+'\t'+mat.split('\t')[2].split(';')[0]+'\t'+str(reference_alleles[mat.split('\t')[2].split(';')[0]][0])+'\t'+str(reference_alleles[mat.split('\t')[2].split(';')[0]][1])+'\t'+mat.split('\t')[6]+'\t'+mat.split('\t')[5]+'\t'+'\t'.join(mat.split('\t')[7:]))
		else:
			out_ase.write(mat)

out_ase.close()

ase_pat=open('ASE.pat.alleles.csv','r')

out_ase=open('ASE.pat.refAllele.csv','w')

header=next(ase_pat)
out_ase.write(header)
for pat in ase_pat:
	if pat.split('\t')[2].split(';')[0] in reference_alleles:
		if pat.split('\t')[3] != reference_alleles[pat.split('\t')[2].split(';')[0]][0]:
			out_ase.write(pat.split('\t')[0]+'\t'+str(reference_alleles[pat.split('\t')[2].split(';')[0]][2])+'\t'+pat.split('\t')[2].split(';')[0]+'\t'+str(reference_alleles[pat.split('\t')[2].split(';')[0]][0])+'\t'+str(reference_alleles[pat.split('\t')[2].split(';')[0]][1])+'\t'+pat.split('\t')[6]+'\t'+pat.split('\t')[5]+'\t'+'\t'.join(pat.split('\t')[7:]))
		else:
			out_ase.write(pat)

out_ase.close()


#################################


ase_mat = read.table('ASE.mat.refAllele.csv',head=T,stringsAsFactors=F)
ase_pat = read.table('ASE.pat.refAllele.csv',head=T,stringsAsFactors=F)

merged_ase=data.frame()

for(i in 1:dim(ase_mat)[1]){

	index=which(ase_pat$variantID == ase_mat$variantID[i])
	if (length(index) < 1){

		merged_ase=rbind(merged_ase,ase_mat[i,])

	} else{

		c_refCount      = ase_mat$refCount[i]+ase_pat$refCount[index]
		c_altCount      = ase_mat$altCount[i]+ase_pat$altCount[index]
		c_totalCount    = ase_mat$totalCount[i]+ase_pat$totalCount[index]
		c_lowMAPQDepth  = ase_mat$lowMAPQDepth[i]+ase_pat$lowMAPQDepth[index]
		c_lowBaseQDepth = ase_mat$lowBaseQDepth[i]+ase_pat$lowBaseQDepth[index]
		c_rawDepth      = ase_mat$rawDepth[i]+ase_pat$rawDepth[index]
		c_otherBases    = ase_mat$otherBases[i]+ase_pat$otherBases[index]
		c_improperPairs = ase_mat$improperPairs[i]+ase_pat$improperPairs[index]
		mod=data.frame(contig=ase_mat$contig[i],position=ase_mat$position[i],variantID=ase_mat$variantID[i],
			       refAllele=ase_mat$refAllele[i],altAllele=ase_mat$altAllele[i],refCount=c_refCount,
			       altCount=c_altCount,totalCount=c_totalCount,lowMAPQDepth=c_lowMAPQDepth,
			       lowBaseQDepth=c_lowBaseQDepth,rawDepth=c_rawDepth,otherBases=c_otherBases,
			       improperPairs=c_improperPairs)
		merged_ase      = rbind(merged_ase,mod)
	}

}

for(i in 1:dim(ase_pat)[1]){

	index=which(ase_mat$variantID == ase_pat$variantID[i])
	if (length(index) < 1){

		merged_ase=rbind(merged_ase,ase_pat[i,])
	}
}


