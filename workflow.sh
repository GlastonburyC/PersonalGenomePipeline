# This should all be excuted from the /data/ directory.
THREAD_NO=8

java -jar ../software/vcf2diploid_v0.2.6a/vcf2diploid.jar -id "$SAMPLE_ID" -chr hg19/hg19.fa -vcf "$SAMPLE_ID"/"$SAMPLE_ID".vcf.gz -outDir $SAMPLE_ID

mv "$SAMPLE_ID"/*_"$SAMPLE_ID"_maternal.fa "$SAMPLE_ID"/maternal/
mv "$SAMPLE_ID"/*_"$SAMPLE_ID"_paternal.fa "$SAMPLE_ID"/paternal/

# Make BAMS sorted by read name
../software/samtools-1.3.1/samtools sort -n "$SAMPLE_ID"/"$SAMPLE_ID"_sorted.bam -o "$SAMPLE_ID"/"$SAMPLE_ID"_sorted.bam.sorted

# convert BAMS to fastq 
../software/bedtools2/bin/bedtools bamtofastq -i "$SAMPLE_ID"/"$SAMPLE_ID"_sorted.bam.sorted -fq "$SAMPLE_ID"/"$SAMPLE_ID".f2.fq -fq2 "$SAMPLE_ID"/"$SAMPLE_ID".f1.fq

# Add check to see whether both ref and personal alignments should be done, or just one (i.e. individuals not in UK10K)
trim_galore -stringency 5 -q 1 -o "$SAMPLE_ID" --phred33 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --paired "$SAMPLE_ID"/"$SAMPLE_ID".f1.fq "$SAMPLE_ID"/"$SAMPLE_ID".f2.fq

# Add trimming and adapter removal step.
perl ../software/prinseq-lite-0.20.4/prinseq-lite.pl -fastq "$SAMPLE_ID"/"$SAMPLE_ID".f1_val_1.fq -fastq2 "$SAMPLE_ID"/"$SAMPLE_ID".f2_val_2.fq -out_good "$SAMPLE_ID"/"$SAMPLE_ID" -trim_tail_left 5 -trim_tail_right 5 -min_len 20

rm "$SAMPLE_ID"/"$SAMPLE_ID".f2_val_2.fq
rm "$SAMPLE_ID"/"$SAMPLE_ID".f1_val_1.fq
rm "$SAMPLE_ID"/"$SAMPLE_ID".f2.fq_trimming_report.txt
rm "$SAMPLE_ID"/"$SAMPLE_ID".f1.fq_trimming_report.txt
rm "$SAMPLE_ID"/"$SAMPLE_ID"_1_singletons.fastq
rm "$SAMPLE_ID"/"$SAMPLE_ID"_2_singletons.fastq
rm "$SAMPLE_ID"/"$SAMPLE_ID"_prinseq_bad*
rm "$SAMPLE_ID"/"$SAMPLE_ID".f1.fq
rm "$SAMPLE_ID"/"$SAMPLE_ID".f2.fq
rm "$SAMPLE_ID"/"$SAMPLE_ID"_sorted.bam.sorted # Ensure this is the name after replacing BAM_IDs with TwinUK IDs.
# Maternal genome generation (suffix arrays etc)
../software/STAR/bin/Linux_x86_64/STAR --runThreadN $THREAD_NO --runMode genomeGenerate --genomeDir "$SAMPLE_ID"/maternal --genomeFastaFiles "$SAMPLE_ID"/maternal/chr*_"$SAMPLE_ID"_maternal.fa

#Paternal genome generation (suffix arrays etc)
../software/STAR/bin/Linux_x86_64/STAR --runThreadN $THREAD_NO --runMode genomeGenerate --genomeDir "$SAMPLE_ID"/paternal --genomeFastaFiles "$SAMPLE_ID"/paternal/chr*_"$SAMPLE_ID"_paternal.fa

#Align paternal
../software/STAR/bin/Linux_x86_64/STAR --runThreadN $THREAD_NO --runMode alignReads --readFilesIn "$SAMPLE_ID"/"$SAMPLE_ID"_1.fastq "$SAMPLE_ID"/"$SAMPLE_ID"_2.fastq --genomeDir "$SAMPLE_ID"/paternal --outSAMstrandField intronMotif --outFilterMultimapNmax 30 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 15 --outMultimapperOrder Random --outSAMunmapped Within --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --sjdbOverhang 48 --outFilterMismatchNmax 6 --outSAMattributes NH nM NM MD HI --outSAMattrRGline  ID:"$SAMPLE_ID"_paternal PU:Illumina PL:Illumina LB:"$SAMPLE_ID"_paternal SM:"$SAMPLE_ID"_paternal CN:Seq_centre --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$SAMPLE_ID"_pat.
mv "$SAMPLE_ID"_pat.Aligned.sortedByCoord.out.bam "$SAMPLE_ID"/paternal/"$SAMPLE_ID"_pat.Aligned.sortedByCoord.out.bam
mv "$SAMPLE_ID"_pat.Chimeric.out.junction "$SAMPLE_ID"/paternal/"$SAMPLE_ID"_pat.Chimeric.out.junction
mv "$SAMPLE_ID"_pat.Chimeric.out.sam "$SAMPLE_ID"/paternal/"$SAMPLE_ID"_pat.Chimeric.out.sam
mv "$SAMPLE_ID"_pat.Log.final.out "$SAMPLE_ID"/paternal/"$SAMPLE_ID"_pat.Log.final.out
mv "$SAMPLE_ID"_pat.SJ.out.tab "$SAMPLE_ID"/paternal/"$SAMPLE_ID"_pat.SJ.out.tab
rm "$SAMPLE_ID"_pat.Log.out
rm "$SAMPLE_ID"_pat.Log.progress.out
#Align maternal
../software/STAR/bin/Linux_x86_64/STAR --runThreadN $THREAD_NO --runMode alignReads --readFilesIn "$SAMPLE_ID"/"$SAMPLE_ID"_1.fastq "$SAMPLE_ID"/"$SAMPLE_ID"_2.fastq --genomeDir "$SAMPLE_ID"/maternal --outSAMstrandField intronMotif --outFilterMultimapNmax 30 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 15 --outMultimapperOrder Random --outSAMunmapped Within --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --sjdbOverhang 48 --outFilterMismatchNmax 6 --outSAMattributes NH nM NM MD HI --outSAMattrRGline  ID:"$SAMPLE_ID"_maternal PU:Illumina PL:Illumina LB:"$SAMPLE_ID"_maternal SM:"$SAMPLE_ID"_maternal CN:Seq_centre --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$SAMPLE_ID"_mat.
mv "$SAMPLE_ID"_mat.Aligned.sortedByCoord.out.bam "$SAMPLE_ID"/maternal/"$SAMPLE_ID"_mat.Aligned.sortedByCoord.out.bam
mv "$SAMPLE_ID"_mat.Chimeric.out.junction "$SAMPLE_ID"/maternal/"$SAMPLE_ID"_mat.Chimeric.out.junction
mv "$SAMPLE_ID"_mat.Chimeric.out.sam "$SAMPLE_ID"/maternal/"$SAMPLE_ID"_mat.Chimeric.out.sam
mv "$SAMPLE_ID"_mat.Log.final.out "$SAMPLE_ID"/maternal/"$SAMPLE_ID"_mat.Log.final.out
mv "$SAMPLE_ID"_mat.SJ.out.tab "$SAMPLE_ID"/maternal/"$SAMPLE_ID"_mat.SJ.out.tab
rm "$SAMPLE_ID"_mat.Log.out
rm "$SAMPLE_ID"_mat.Log.progress.out
# Add standard reference alignment too.

../software/STAR/bin/Linux_x86_64/STAR --runThreadN $THREAD_NO --runMode alignReads --readFilesIn "$SAMPLE_ID"/"$SAMPLE_ID"_1.fastq "$SAMPLE_ID"/"$SAMPLE_ID"_2.fastq --genomeDir hg19 --outSAMstrandField intronMotif --outFilterMultimapNmax 30 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 15 --outMultimapperOrder Random --outSAMunmapped Within --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --sjdbOverhang 48 --outFilterMismatchNmax 6 --outSAMattributes NH nM NM MD HI --outSAMattrRGline  ID:"$SAMPLE_ID"_maternal PU:Illumina PL:Illumina LB:"$SAMPLE_ID"_maternal SM:"$SAMPLE_ID"_maternal CN:Seq_centre --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$SAMPLE_ID"_ref.
mv "$SAMPLE_ID"_ref.Aligned.sortedByCoord.out.bam "$SAMPLE_ID"/reference/"$SAMPLE_ID"_ref.Aligned.sortedByCoord.out.bam
mv  "$SAMPLE_ID"_ref.Chimeric.out.junction "$SAMPLE_ID"/reference/"$SAMPLE_ID"_ref.Chimeric.out.junction
mv  "$SAMPLE_ID"_ref.Chimeric.out.sam "$SAMPLE_ID"/reference/"$SAMPLE_ID"_ref.Chimeric.out.sam
mv  "$SAMPLE_ID"_ref.Log.final.out "$SAMPLE_ID"/reference/"$SAMPLE_ID"_ref.Log.final.out
mv  "$SAMPLE_ID"_ref.SJ.out.tab "$SAMPLE_ID"/reference/"$SAMPLE_ID"_ref.SJ.out.tab
rm "$SAMPLE_ID"_ref.Log.out
rm "$SAMPLE_ID"_ref.Log.progress.out
########################################
#
## Produce consensus bams, in which the best read per haplotype is selected.
python ../software/PersonalGenomePipeline/pipeline.seperateBAMs.py 1152

# Sort the VCFs else picard/GATK throws a fit.
../software/vcftools/src/perl/vcf-sort -c "$SAMPLE_ID"/"$SAMPLE_ID".hets.phased.vcf.gz > "$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.vcf

pigz --best -k "$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.vcf

rm "$SAMPLE_ID"/"$SAMPLE_ID".hets.phased.vcf.gz
rm "$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.vcf
#########################################
# run BASH script to concatenate reference genomes, and rename chain 1_maternal > 1

cp "$SAMPLE_ID"/paternal.chain "$SAMPLE_ID"/"$SAMPLE_ID".paternal.edit.chain
cp "$SAMPLE_ID"/maternal.chain "$SAMPLE_ID"/"$SAMPLE_ID".maternal.edit.chain
for i in {1..22}
do
	sed -i 's/'"$i"'_maternal/'"$i"'/g' "$SAMPLE_ID"/"$SAMPLE_ID".maternal.edit.chain
	sed -i 's/'"$i"'_paternal/'"$i"'/g' "$SAMPLE_ID"/"$SAMPLE_ID".paternal.edit.chain
	cat "$SAMPLE_ID"/paternal/"$i"_"$SAMPLE_ID"_paternal.fa >> "$SAMPLE_ID"/paternal/"$SAMPLE_ID".all_paternal.fa
	cat "$SAMPLE_ID"/maternal/"$i"_"$SAMPLE_ID"_maternal.fa >> "$SAMPLE_ID"/paternal/"$SAMPLE_ID".all_maternal.fa

done

sed -i 's/X_paternal/X/g' "$SAMPLE_ID"/"$SAMPLE_ID".paternal.edit.chain
sed -i 's/Y_paternal/Y/g' "$SAMPLE_ID"/"$SAMPLE_ID".paternal.edit.chain
sed -i 's/MT_paternal/MT/g' "$SAMPLE_ID"/"$SAMPLE_ID".paternal.edit.chain

sed -i 's/X_maternal/X/g' "$SAMPLE_ID"/"$SAMPLE_ID".maternal.edit.chain
sed -i 's/Y_maternal/Y/g' "$SAMPLE_ID"/"$SAMPLE_ID".maternal.edit.chain
sed -i 's/MT_maternal/MT/g' "$SAMPLE_ID"/"$SAMPLE_ID".maternal.edit.chain

cat X_131_paternal.fa >> "$SAMPLE_ID"/"$SAMPLE_ID".all_paternal.fa
cat Y_131_paternal.fa >> "$SAMPLE_ID"/"$SAMPLE_ID".all_paternal.fa
cat MT_131_paternal.fa >> "$SAMPLE_ID"/"$SAMPLE_ID".all_paternal.fa

cat X_131_maternal.fa >> "$SAMPLE_ID"/"$SAMPLE_ID".all_maternal.fa
cat Y_131_maternal.fa >> "$SAMPLE_ID"/"$SAMPLE_ID".all_maternal.fa
cat MT_131_maternal.fa >> "$SAMPLE_ID"/"$SAMPLE_ID".all_maternal.fa

#################################

# rename maternal and paternal fasta files with correct scaffold/chr name.
python renameFaChr.py "$SAMPLE_ID"/"$SAMPLE_ID".all_maternal.fa "$SAMPLE_ID"/"$SAMPLE_ID".maternal.renamed.fa "$SAMPLE_ID"/"$SAMPLE_ID".all_maternal.fa "$SAMPLE_ID"/"$SAMPLE_ID".paternal.renamed.fa

# Remap VCF to parental genome coordinates using modified CrossMap (fixedBugs)
python mod.py vcf "$SAMPLE_ID"/"$SAMPLE_ID".maternal.edit.chain "$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.renamed.vcf "$SAMPLE_ID"/"$SAMPLE_ID".maternal.renamed.fa "$SAMPLE_ID"/"$SAMPLE_ID".maternal.vcf

python mod.py vcf "$SAMPLE_ID"/"$SAMPLE_ID".paternal.edit.chain "$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.renamed.vcf "$SAMPLE_ID"/"$SAMPLE_ID".paternal.renamed.fa "$SAMPLE_ID"/"$SAMPLE_ID".paternal.vcf

# This script swaps the REF and ALT alleles according to whether it's the maternal or paternal haplotype.
haplotypeVCFAlleles.py "$SAMPLE_ID"/"$SAMPLE_ID".maternal.vcf "$SAMPLE_ID"/"$SAMPLE_ID".paternal.vcf

/usr/lib/jvm/java-8-oracle/bin/java -jar picard-tools-2.1.1/picard.jar CreateSequenceDictionary R=""$SAMPLE_ID"/$SAMPLE_ID".paternal.renamed.fa O="$SAMPLE_ID"/"$SAMPLE_ID".paternal.renamed.dict
/usr/lib/jvm/java-8-oracle/bin/java -jar picard-tools-2.1.1/picard.jar CreateSequenceDictionary R=""$SAMPLE_ID"/$SAMPLE_ID".maternal.renamed.fa O=""$SAMPLE_ID"/$SAMPLE_ID".maternal.renamed.dict

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
   and for easier interpretation, its better to have ASE in terms of universal REF / ALT. i.e. adding up all the REF/ALT counts 
   for each haplotype, consolidated into a single ASE file per individual.'''

python ASERefCord.py "$SAMPLE_ID".ASE.mat.csv "$SAMPLE_ID".ASE.pat.csv "$SAMPLE_ID".hets.GATK.sorted.vcf "$SAMPLE_ID".maternal.vcf "$SAMPLE_ID".paternal.vcf "$SAMPLE_ID".maternal.alleles.csv "$SAMPLE_ID".paternal.alleles.csv "$SAMPLE_ID".maternal.ref.csv "$SAMPLE_ID".paternal.ref.csv

#################################

R CMD BATCH ASERefOut.R --args "$SAMPLE_ID".maternal.ref.csv "$SAMPLE_ID".paternal.ref.csv "$SAMPLE_ID".ASE.csv
