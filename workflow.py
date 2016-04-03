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
	#sed -i 's/'"$i"'_maternal/'"$i"'/g' maternal.edit.chain
	#sed -i 's/'"$i"'_paternal/'"$i"'/g' paternal.edit.chain
	#cat "$i"_131_paternal.fa >> all_paternal.fa
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

fasta= open('all_paternal.fa')
newnames= open('rename.txt')
newfasta= open('paternal.renamed.fa', 'w')

for line in fasta:
    if line.startswith('>'):
        newname= newnames.readline()
        newfasta.write(newname)
    else:
        newfasta.write(line)

fasta.close()
newnames.close()
newfasta.close()

fasta= open('all_maternal.fa')
newnames= open('rename.txt')
newfasta= open('maternal.renamed.fa', 'w')

for line in fasta:
    if line.startswith('>'):
        newname= newnames.readline()
        newfasta.write(newname)
    else:
        newfasta.write(line)

fasta.close()
newnames.close()
newfasta.close()

# Remap VCF to parental genome coordinates using modified CrossMap (fixedBugs)
sudo python mod.py vcf maternal.edit.chain 131.hets.GATK.sorted.vcf maternal.renamed.fa 131.maternal.vcf

sudo python mod.py vcf paternal.edit.chain 131.hets.GATK.sorted.vcf paternal.renamed.fa 131.paternal.vcf






