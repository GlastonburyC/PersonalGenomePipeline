module load samtools
module load vcftools
module load R  # change depending on your cluster environment
module load python/2.7
# Maternal genome generation (suffix arrays etc)
STAR/bin/Linux_x86_64/STAR --runThreadN 8 --runMode genomeGenerate --genomeDir "$SAMPLE_ID"/maternal/ --genomeFastaFiles *_"$SAMPLE_ID"_maternal.fa

#Paternal genome generation (suffix arrays etc)
STAR/bin/Linux_x86_64/STAR --runThreadN 8 --runMode genomeGenerate --genomeDir "$SAMPLE_ID"/paternal/ --genomeFastaFiles *_"$SAMPLE_ID"_paternal.fa

#Align paternal
STAR/bin/Linux_x86_64/STAR --runThreadN 8 --runMode alignReads --readFilesIn "$SAMPLE_ID"/F226_new_120618_7_sorted.bam.f1_val_1.fq "$SAMPLE_ID"/F226_new_120618_7_sorted.bam.f2_val_2.fq --genomeDir "$SAMPLE_ID"/maternal --outFilterMultimapNmax 30 --sjdbOverhang 48 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattrIHstart 0 --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --chimSegmentMin 15 --outMultimapperOrder Random --outFilterMismatchNmax 4 --outSAMunmapped Within --outSAMattributes NH HI NM MD

#Align maternal
STAR/bin/Linux_x86_64/STAR --runThreadN 8 --runMode alignReads --readFilesIn "$SAMPLE_ID"/F226_new_120618_7_sorted.bam.f1_val_1.fq "$SAMPLE_ID"/F226_new_120618_7_sorted.bam.f2_val_2.fq --genomeDir "$SAMPLE_ID"/maternal --outFilterMultimapNmax 30 --sjdbOverhang 48 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattrIHstart 0 --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --chimSegmentMin 15 --outMultimapperOrder Random --outFilterMismatchNmax 4 --outSAMunmapped Within --outSAMattributes NH HI NM MD

########################################

vcf-sort -c "$SAMPLE_ID".hets.GATK.vcf > "$SAMPLE_ID".hets.GATK.sorted.vcf

#########################################
# run BASH script to concatenate reference genomes, and rename chain 1_maternal > 1

cp "$SAMPLE_ID".paternal.chain "$SAMPLE_ID".paternal.edit.chain
cp "$SAMPLE_ID".maternal.chain "$SAMPLE_ID".maternal.edit.chain
for i in {1..22}
do
	sed -i 's/'"$i"'_maternal/'"$i"'/g' "$SAMPLE_ID".maternal.edit.chain
	sed -i 's/'"$i"'_paternal/'"$i"'/g' "$SAMPLE_ID".paternal.edit.chain
	cat "$i"_"$SAMPLE_ID"_paternal.fa >> "$SAMPLE_ID".all_paternal.fa
	cat "$i"_"$SAMPLE_ID"_maternal.fa >> "$SAMPLE_ID".all_maternal.fa

done

sed -i 's/X_paternal/X/g' "$SAMPLE_ID".paternal.edit.chain
sed -i 's/Y_paternal/Y/g' "$SAMPLE_ID".paternal.edit.chain
sed -i 's/MT_paternal/MT/g' "$SAMPLE_ID".paternal.edit.chain

sed -i 's/X_maternal/X/g' "$SAMPLE_ID".maternal.edit.chain
sed -i 's/Y_maternal/Y/g' "$SAMPLE_ID".maternal.edit.chain
sed -i 's/MT_maternal/MT/g' "$SAMPLE_ID".maternal.edit.chain

cat X_131_paternal.fa >> "$SAMPLE_ID".all_paternal.fa
cat Y_131_paternal.fa >> "$SAMPLE_ID".all_paternal.fa
cat MT_131_paternal.fa >> "$SAMPLE_ID".all_paternal.fa

cat X_131_maternal.fa >> "$SAMPLE_ID".all_maternal.fa
cat Y_131_maternal.fa >> "$SAMPLE_ID".all_maternal.fa
cat MT_131_maternal.fa >> "$SAMPLE_ID".all_maternal.fa

#################################

# rename maternal and paternal fasta files with correct scaffold/chr name.
python renameFaChr.py "$SAMPLE_ID".all_maternal.fa "$SAMPLE_ID".maternal.renamed.fa "$SAMPLE_ID".all_maternal.fa "$SAMPLE_ID".paternal.renamed.fa

# Remap VCF to parental genome coordinates using modified CrossMap (fixedBugs)
python mod.py vcf "$SAMPLE_ID".maternal.edit.chain "$SAMPLE_ID".hets.GATK.sorted.renamed.vcf "$SAMPLE_ID".maternal.renamed.fa "$SAMPLE_ID".maternal.vcf

python mod.py vcf "$SAMPLE_ID".paternal.edit.chain "$SAMPLE_ID".hets.GATK.sorted.renamed.vcf "$SAMPLE_ID".paternal.renamed.fa "$SAMPLE_ID".paternal.vcf

# This script swaps the REF and ALT alleles according to whether it's the maternal or paternal haplotype.
haplotypeVCFAlleles.py "$SAMPLE_ID".maternal.vcf "$SAMPLE_ID".paternal.vcf

/usr/lib/jvm/java-8-oracle/bin/java -jar picard-tools-2.1.1/picard.jar CreateSequenceDictionary R="$SAMPLE_ID".paternal.renamed.fa O="$SAMPLE_ID".paternal.renamed.dict
/usr/lib/jvm/java-8-oracle/bin/java -jar picard-tools-2.1.1/picard.jar CreateSequenceDictionary R="$SAMPLE_ID".maternal.renamed.fa O="$SAMPLE_ID".maternal.renamed.dict

/usr/lib/jvm/java-8-oracle/bin/java -jar picard-tools-2.1.1/picard.jar AddOrReplaceReadGroups I="$SAMPLE_ID"/consensus.mat.filtered.sorted.bam O="$SAMPLE_ID"/consensus.mat.filtered.sorted.readGroup.bam  RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
/usr/lib/jvm/java-8-oracle/bin/java -jar picard-tools-2.1.1/picard.jar AddOrReplaceReadGroups I="$SAMPLE_ID"/consensus.pat.filtered.sorted.bam O="$SAMPLE_ID"/consensus.pat.filtered.sorted.readGroup.bam  RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20

# both say maternal for chromosome because maternal was used as a template BAM when selecting best reads.
samtools view -h "$SAMPLE_ID"/consensus.mat.filtered.sorted.readGroup.bam | sed -e 's/_maternal//g' >> "$SAMPLE_ID"/consensus.mat.filtered.sorted2.bam
samtools view -h "$SAMPLE_ID"/consensus.pat.filtered.sorted.readGroup.bam | sed -e 's/_maternal//g' >> "$SAMPLE_ID"/consensus.pat.filtered.sorted2.bam

mv "$SAMPLE_ID"/consensus.mat.filtered.sorted2.bam "$SAMPLE_ID"/consensus.mat.filtered.sorted.readGroup.bam
mv "$SAMPLE_ID"/consensus.pat.filtered.sorted2.bam "$SAMPLE_ID"/consensus.pat.filtered.sorted.readGroup.bam

/usr/lib/jvm/java-8-oracle/bin/java -jar picard-tools-2.1.1/picard.jar ReorderSam I="$SAMPLE_ID"/consensus.mat.filtered.sorted.readGroup.bam O="$SAMPLE_ID"/consensus.mat.filtered.sorted.readGroup.sorted.bam R="$SAMPLE_ID".maternal.renamed.fa ALLOW_CONTIG_LENGTH_DISCORDANCE=true
/usr/lib/jvm/java-8-oracle/bin/java -jar picard-tools-2.1.1/picard.jar ReorderSam I="$SAMPLE_ID"/consensus.pat.filtered.sorted.readGroup.bam O="$SAMPLE_ID"/consensus.pat.filtered.sorted.readGroup.sorted.bam R="$SAMPLE_ID".paternal.renamed.fa ALLOW_CONTIG_LENGTH_DISCORDANCE=true

java -jar GenomeAnalysisTK.jar -R "$SAMPLE_ID".maternal.renamed.fa -T ASEReadCounter -o "$SAMPLE_ID".ASE.mat.csv -I "$SAMPLE_ID"/consensus.mat.filtered.sorted.readGroup.sorted.bam -sites "$SAMPLE_ID".maternal.vcf -dels -U ALLOW_N_CIGAR_READS -S SILENT
java -jar GenomeAnalysisTK.jar -R "$SAMPLE_ID".paternal.renamed.fa -T ASEReadCounter -o "$SAMPLE_ID".ASE.pat.csv -I "$SAMPLE_ID"/consensus.pat.filtered.sorted.readGroup.sorted.bam "$SAMPLE_ID".paternal.vcf -dels -U ALLOW_N_CIGAR_READS -S SILENT

'''Currently, the ASE output is relative to the haplotype reference used (Maternal reference or paternal reference) 
   i.e. REF in ASE.mat.csv = maternal allele etc.To assess REF Bias (which we should not expect due to using Personal genomes), 
   and for easier interpretation, its better to have ASE in terms of universal REF / ALT. i.e. all adding up all the REF/ALT counts 
   for each haplotype, consolidated into a single ASE file per individual.'''


python ASERefCord.py "$SAMPLE_ID".ASE.mat.csv "$SAMPLE_ID".ASE.pat.csv "$SAMPLE_ID".hets.GATK.sorted.vcf "$SAMPLE_ID".maternal.vcf "$SAMPLE_ID".paternal.vcf "$SAMPLE_ID".maternal.alleles.csv "$SAMPLE_ID".paternal.alleles.csv "$SAMPLE_ID".maternal.ref.csv "$SAMPLE_ID".paternal.ref.csv

#################################

R CMD BATCH ASERefOut.R --args "$SAMPLE_ID".maternal.ref.csv "$SAMPLE_ID".paternal.ref.csv "$SAMPLE_ID".ASE.csv
