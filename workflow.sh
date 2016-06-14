# USAGE:   sh SCRIPT.sh TEST_FILE CHECK_FILE
# This should all be excuted from the /data/ directory.
THREAD_NO=$3
ALL_SAMPLES=$1
UK10K_SAMPLES=$2

## for each line in TEST_FILE
while read line ; do

#!/bin/bash 
# 
#SBATCH -n 8 
## check if line exist in CHECK_FILE; then assign result to variable
SAMPLE_ID=$(grep "^${line}$" ${UK10K_SAMPLES})


## if variable is blank (meaning ALL_SAMPLES line not found in UK10K)
## print 'false' and exit
if [[ -z $SAMPLE_ID ]] ; then
 
echo '
#!/bin/bash 
# 
#SBATCH -N 1 
# number of nodes 
#SBATCH -n 8 

echo "1. Sorting BAM file by read name."
../software/samtools-1.3.1/samtools sort -n "$SAMPLE_ID"/"$SAMPLE_ID"_sorted.bam -o "$SAMPLE_ID"/"$SAMPLE_ID"_sorted.bam.sorted
echo "2. BAM > fastq."
../software/bedtools2/bin/bedtools bamtofastq -i "$SAMPLE_ID"/"$SAMPLE_ID"_sorted.bam.sorted -fq "$SAMPLE_ID"/"$SAMPLE_ID".f2.fq -fq2 "$SAMPLE_ID"/"$SAMPLE_ID".f1.fq
echo "3. Trimming adapters."
../software/trim_galore_zip/trim_galore -stringency 5 -q 1 -o "$SAMPLE_ID" --phred33 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --paired "$SAMPLE_ID"/"$SAMPLE_ID".f1.fq "$SAMPLE_ID"/"$SAMPLE_ID".f2.fq
echo "4. Trimming polyA tails."
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
rm "$SAMPLE_ID"/"$SAMPLE_ID"_sorted.bam.sorted
echo "5. Aligning transcriptome with STAR."
../software/STAR/bin/Linux_x86_64/STAR --runThreadN $THREAD_NO --runMode alignReads --readFilesIn "$SAMPLE_ID"/"$SAMPLE_ID"_1.fastq "$SAMPLE_ID"/"$SAMPLE_ID"_2.fastq --genomeDir hg19 --outSAMstrandField intronMotif --outFilterMultimapNmax 30 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 15 --outMultimapperOrder Random --outSAMunmapped Within --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --sjdbOverhang 48 --outFilterMismatchNmax 6 --outSAMattributes NH nM NM MD HI --outSAMattrRGline  ID:"$SAMPLE_ID"_maternal PU:Illumina PL:Illumina LB:"$SAMPLE_ID"_maternal SM:"$SAMPLE_ID"_maternal CN:Seq_centre --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$SAMPLE_ID"_ref.

rm -r "$SAMPLE_ID"_mat._STARtmp
mv "$SAMPLE_ID"_ref.Aligned.sortedByCoord.out.bam "$SAMPLE_ID"/reference/"$SAMPLE_ID"_ref.Aligned.sortedByCoord.out.bam
mv "$SAMPLE_ID"_ref.Chimeric.out.junction "$SAMPLE_ID"/reference/"$SAMPLE_ID"_ref.Chimeric.out.junction
mv "$SAMPLE_ID"_ref.Chimeric.out.sam "$SAMPLE_ID"/reference/"$SAMPLE_ID"_ref.Chimeric.out.sam
mv "$SAMPLE_ID"_ref.Log.final.out "$SAMPLE_ID"/reference/"$SAMPLE_ID"_ref.Log.final.out
mv "$SAMPLE_ID"_ref.SJ.out.tab "$SAMPLE_ID"/reference/"$SAMPLE_ID"_ref.SJ.out.tab
rm "$SAMPLE_ID"_ref.Log.out
rm "$SAMPLE_ID"_ref.Log.progress.out
echo "5. Filtering aligned BAM QUAL > 30, properly paired"
../software/samtools-1.3.1/samtools view -b -F4 -q 30 "$SAMPLE_ID"/reference/"$SAMPLE_ID"_ref.Aligned.sortedByCoord.out.bam -o "$SAMPLE_ID"/reference/"$SAMPLE_ID".filtered.bam
../software/subread-1.5.0-p3-Linux-x86_64/bin/featureCounts -p -T 8 -a gencode.v19.annotation.gtf -o "$SAMPLE_ID"/reference/"$SAMPLE_ID".GeneCount_Ref.txt "$SAMPLE_ID"/reference/"$SAMPLE_ID".filtered.bam

echo "6. Finished Successfully"
rm "$SAMPLE_ID"/"$SAMPLE_ID"_sorted.bam.sorted
rm "$SAMPLE_ID"/"$SAMPLE_ID"_sorted.bam ' > $SAMPLE_ID.refOnly.sh
else
echo "
#!/bin/bash 
# 
#SBATCH -N 1 
# number of nodes 
#SBATCH -n 8 

echo '1. Editing reference genome (hg19) with maternal and paternal variants (phased, QUAL > 30).'
java -jar ../software/vcf2diploid_v0.2.6a/vcf2diploid.jar -id "$SAMPLE_ID" -chr hg19/hg19.fa -vcf "$SAMPLE_ID"/"$SAMPLE_ID".vcf.gz -outDir $SAMPLE_ID

mv "$SAMPLE_ID"/*_"$SAMPLE_ID"_maternal.fa "$SAMPLE_ID"/maternal/
mv "$SAMPLE_ID"/*_"$SAMPLE_ID"_paternal.fa "$SAMPLE_ID"/paternal/

# Make BAMS sorted by read name
echo '2. Sorting BAM by readname.'
../software/samtools-1.3.1/samtools sort -n "$SAMPLE_ID"/"$SAMPLE_ID"_sorted.bam -o "$SAMPLE_ID"/"$SAMPLE_ID"_sorted.bam.sorted

# convert BAMS to fastq 
echo '3. BAM > fastq.'
../software/bedtools2/bin/bedtools bamtofastq -i "$SAMPLE_ID"/"$SAMPLE_ID"_sorted.bam.sorted -fq "$SAMPLE_ID"/"$SAMPLE_ID".f2.fq -fq2 "$SAMPLE_ID"/"$SAMPLE_ID".f1.fq

# Add check to see whether both ref and personal alignments should be done, or just one (i.e. individuals not in UK10K)
echo '4. Trimming adapters.'
../software/trim_galore_zip/trim_galore -stringency 5 -q 1 -o "$SAMPLE_ID" --phred33 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --paired "$SAMPLE_ID"/"$SAMPLE_ID".f1.fq "$SAMPLE_ID"/"$SAMPLE_ID".f2.fq
# Add trimming and adapter removal step.
echo '5. Trimming polyA tails.'
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
echo '5. Generating Maternal reference suffix array.'
../software/STAR/bin/Linux_x86_64/STAR --runThreadN $THREAD_NO --runMode genomeGenerate --genomeDir "$SAMPLE_ID"/maternal --genomeFastaFiles "$SAMPLE_ID"/maternal/chr*_"$SAMPLE_ID"_maternal.fa

#Paternal genome generation (suffix arrays etc)
echo '6. Generating Paternal reference suffix array.'
../software/STAR/bin/Linux_x86_64/STAR --runThreadN $THREAD_NO --runMode genomeGenerate --genomeDir "$SAMPLE_ID"/paternal --genomeFastaFiles "$SAMPLE_ID"/paternal/chr*_"$SAMPLE_ID"_paternal.fa

#Align paternal
echo '7. Aligning to paternally editted reference genome (STAR).'
../software/STAR/bin/Linux_x86_64/STAR --runThreadN $THREAD_NO --runMode alignReads --readFilesIn "$SAMPLE_ID"/"$SAMPLE_ID"_1.fastq "$SAMPLE_ID"/"$SAMPLE_ID"_2.fastq --genomeDir "$SAMPLE_ID"/paternal --outSAMstrandField intronMotif --outFilterMultimapNmax 30 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 15 --outMultimapperOrder Random --outSAMunmapped Within --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --sjdbOverhang 48 --outFilterMismatchNmax 6 --outSAMattributes NH nM NM MD HI --outSAMattrRGline  ID:"$SAMPLE_ID"_paternal PU:Illumina PL:Illumina LB:"$SAMPLE_ID"_paternal SM:"$SAMPLE_ID"_paternal CN:Seq_centre --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$SAMPLE_ID"_pat.
mv "$SAMPLE_ID"_pat.Aligned.sortedByCoord.out.bam "$SAMPLE_ID"/paternal/"$SAMPLE_ID"_pat.Aligned.sortedByCoord.out.bam
mv "$SAMPLE_ID"_pat.Chimeric.out.junction "$SAMPLE_ID"/paternal/"$SAMPLE_ID"_pat.Chimeric.out.junction
mv "$SAMPLE_ID"_pat.Chimeric.out.sam "$SAMPLE_ID"/paternal/"$SAMPLE_ID"_pat.Chimeric.out.sam
mv "$SAMPLE_ID"_pat.Log.final.out "$SAMPLE_ID"/paternal/"$SAMPLE_ID"_pat.Log.final.out
mv "$SAMPLE_ID"_pat.SJ.out.tab "$SAMPLE_ID"/paternal/"$SAMPLE_ID"_pat.SJ.out.tab
rm "$SAMPLE_ID"_pat.Log.out
rm "$SAMPLE_ID"_pat.Log.progress.out
#Align maternal
echo '8. Aligning to maternally editted reference genome (STAR).'
../software/STAR/bin/Linux_x86_64/STAR --runThreadN $THREAD_NO --runMode alignReads --readFilesIn "$SAMPLE_ID"/"$SAMPLE_ID"_1.fastq "$SAMPLE_ID"/"$SAMPLE_ID"_2.fastq --genomeDir "$SAMPLE_ID"/maternal --outSAMstrandField intronMotif --outFilterMultimapNmax 30 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 15 --outMultimapperOrder Random --outSAMunmapped Within --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --sjdbOverhang 48 --outFilterMismatchNmax 6 --outSAMattributes NH nM NM MD HI --outSAMattrRGline  ID:"$SAMPLE_ID"_maternal PU:Illumina PL:Illumina LB:"$SAMPLE_ID"_maternal SM:"$SAMPLE_ID"_maternal CN:Seq_centre --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$SAMPLE_ID"_mat.
mv "$SAMPLE_ID"_mat.Aligned.sortedByCoord.out.bam "$SAMPLE_ID"/maternal/"$SAMPLE_ID"_mat.Aligned.sortedByCoord.out.bam
mv "$SAMPLE_ID"_mat.Chimeric.out.junction "$SAMPLE_ID"/maternal/"$SAMPLE_ID"_mat.Chimeric.out.junction
mv "$SAMPLE_ID"_mat.Chimeric.out.sam "$SAMPLE_ID"/maternal/"$SAMPLE_ID"_mat.Chimeric.out.sam
mv "$SAMPLE_ID"_mat.Log.final.out "$SAMPLE_ID"/maternal/"$SAMPLE_ID"_mat.Log.final.out
mv "$SAMPLE_ID"_mat.SJ.out.tab "$SAMPLE_ID"/maternal/"$SAMPLE_ID"_mat.SJ.out.tab
rm "$SAMPLE_ID"_mat.Log.out
rm "$SAMPLE_ID"_mat.Log.progress.out

# Add standard reference alignment too.
echo '9. Aligning to hg19 reference genome (STAR).'
../software/STAR/bin/Linux_x86_64/STAR --runThreadN $THREAD_NO --runMode alignReads --readFilesIn "$SAMPLE_ID"/"$SAMPLE_ID"_1.fastq "$SAMPLE_ID"/"$SAMPLE_ID"_2.fastq --genomeDir hg19 --outSAMstrandField intronMotif --outFilterMultimapNmax 30 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 15 --outMultimapperOrder Random --outSAMunmapped Within --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --sjdbOverhang 48 --outFilterMismatchNmax 6 --outSAMattributes NH nM NM MD HI --outSAMattrRGline  ID:"$SAMPLE_ID"_maternal PU:Illumina PL:Illumina LB:"$SAMPLE_ID"_maternal SM:"$SAMPLE_ID"_maternal CN:Seq_centre --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$SAMPLE_ID"_ref.

mv "$SAMPLE_ID"_ref.Aligned.sortedByCoord.out.bam "$SAMPLE_ID"/reference/"$SAMPLE_ID"_ref.Aligned.sortedByCoord.out.bam
mv "$SAMPLE_ID"_ref.Chimeric.out.junction "$SAMPLE_ID"/reference/"$SAMPLE_ID"_ref.Chimeric.out.junction
mv "$SAMPLE_ID"_ref.Chimeric.out.sam "$SAMPLE_ID"/reference/"$SAMPLE_ID"_ref.Chimeric.out.sam
mv "$SAMPLE_ID"_ref.Log.final.out "$SAMPLE_ID"/reference/"$SAMPLE_ID"_ref.Log.final.out
mv "$SAMPLE_ID"_ref.SJ.out.tab "$SAMPLE_ID"/reference/"$SAMPLE_ID"_ref.SJ.out.tab
rm "$SAMPLE_ID"_ref.Log.out
rm "$SAMPLE_ID"_ref.Log.progress.out
########################################
#

rm -r "$SAMPLE_ID"_mat._STARtmp
rm -r "$SAMPLE_ID"_pat._STARtmp
rm -r "$SAMPLE_ID"_ref._STARtmp


echo '10. Filtering alignment (QUAL > 30, properly paired).'
../software/samtools-1.3.1/samtools view -b -F4 -q 30 "$SAMPLE_ID"/reference/"$SAMPLE_ID"_ref.Aligned.sortedByCoord.out.bam -o "$SAMPLE_ID"/reference/"$SAMPLE_ID".filtered.bam

## Produce consensus bams, in which the best read per haplotype is selected.
echo '11. Running seperateBAMs.py - Selects best read alignment from each haplotype based on edit distance'
python ../software/PersonalGenomePipeline/pipeline.seperateBAMs.py "$SAMPLE_ID"
#
#for all variants without an rsid - assign them chrpos.
echo '12. replacing empty rsid tags in VCF with chr-pos.'
python ../software/PersonalGenomePipeline/nameVariants.py "$SAMPLE_ID" "$SAMPLE_ID".hets.phased.vcf.gz "$SAMPLE_ID".hets.phased.vcf.gz2

mv "$SAMPLE_ID"/"$SAMPLE_ID".hets.phased.vcf.gz2 "$SAMPLE_ID"/"$SAMPLE_ID".hets.phased.vcf.gz

# Sort the VCFs else picard/GATK throws a fit.
echo '13. Sorting VCF by position.'
../software/vcftools/src/perl/vcf-sort -c "$SAMPLE_ID"/"$SAMPLE_ID".hets.phased.vcf.gz > "$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.vcf

pigz --best -k "$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.vcf

rm "$SAMPLE_ID"/"$SAMPLE_ID".hets.phased.vcf.gz
rm "$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.vcf
#########################################
# run BASH script to concatenate reference genomes, and rename chain 1_maternal > 1

echo '14. Harmonizing chain files, fasta files and VCF on how they refer to chromosomes. (Step 1).'
cp "$SAMPLE_ID"/paternal.chain "$SAMPLE_ID"/"$SAMPLE_ID".paternal.edit.chain
cp "$SAMPLE_ID"/maternal.chain "$SAMPLE_ID"/"$SAMPLE_ID".maternal.edit.chain
for i in {1..22}
do
	sed -i 's/\$i_maternal/\$i/g' "$SAMPLE_ID"/"$SAMPLE_ID".maternal.edit.chain
	sed -i 's/\$i_paternal/\$i/g' "$SAMPLE_ID"/"$SAMPLE_ID".paternal.edit.chain
	cat "$SAMPLE_ID"/paternal/chr\$i_"$SAMPLE_ID"_paternal.fa >> "$SAMPLE_ID"/paternal/"$SAMPLE_ID".all_paternal.fa
	cat "$SAMPLE_ID"/maternal/chr\$i_"$SAMPLE_ID"_maternal.fa >> "$SAMPLE_ID"/maternal/"$SAMPLE_ID".all_maternal.fa

done
echo '15. Harmonizing chain files, fasta files and VCF on how they refer to chromosomes. (Step 2).'

sed -i 's/X_paternal/X/g' "$SAMPLE_ID"/"$SAMPLE_ID".paternal.edit.chain
#sed -i 's/Y_paternal/Y/g' "$SAMPLE_ID"/"$SAMPLE_ID".paternal.edit.chain
sed -i 's/M_paternal/MT/g' "$SAMPLE_ID"/"$SAMPLE_ID".paternal.edit.chain

sed -i 's/X_maternal/X/g' "$SAMPLE_ID"/"$SAMPLE_ID".maternal.edit.chain
#sed -i 's/Y_maternal/Y/g' "$SAMPLE_ID"/"$SAMPLE_ID".maternal.edit.chain
sed -i 's/M_maternal/MT/g' "$SAMPLE_ID"/"$SAMPLE_ID".maternal.edit.chain

echo '16. Harmonizing chain files, fasta files and VCF on how they refer to chromosomes. (Step 3).'
cat "$SAMPLE_ID"/paternal/chrX_"$SAMPLE_ID"_paternal.fa >> "$SAMPLE_ID"/paternal/"$SAMPLE_ID".all_paternal.fa
#cat "$SAMPLE_ID"/paternal/chrY_"$SAMPLE_ID"_paternal.fa >> "$SAMPLE_ID"/paternal/"$SAMPLE_ID".all_paternal.fa
cat "$SAMPLE_ID"/paternal/chrM_"$SAMPLE_ID"_paternal.fa >> "$SAMPLE_ID"/paternal/"$SAMPLE_ID".all_paternal.fa

cat "$SAMPLE_ID"/maternal/chrX_"$SAMPLE_ID"_maternal.fa >> "$SAMPLE_ID"/maternal/"$SAMPLE_ID".all_maternal.fa
#cat "$SAMPLE_ID"/maternal/chrY_"$SAMPLE_ID"_maternal.fa >> "$SAMPLE_ID"/maternal/"$SAMPLE_ID".all_maternal.fa
cat "$SAMPLE_ID"/maternal/chrM_"$SAMPLE_ID"_maternal.fa >> "$SAMPLE_ID"/maternal/"$SAMPLE_ID".all_maternal.fa

#################################

# rename maternal and paternal fasta files with correct scaffold/chr name.
echo '17. Harmonizing chain files, fasta files and VCF on how they refer to chromosomes. (Step 4).'
python ../software/PersonalGenomePipeline/renameFaChr.py "$SAMPLE_ID"/maternal/"$SAMPLE_ID".all_maternal.fa "$SAMPLE_ID"/maternal/"$SAMPLE_ID".maternal.renamed.fa "$SAMPLE_ID"/paternal/"$SAMPLE_ID".all_paternal.fa "$SAMPLE_ID"/paternal/"$SAMPLE_ID".paternal.renamed.fa

pigz -d "$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.vcf.gz

echo '18. Harmonizing chain files, fasta files and VCF on how they refer to chromosomes. (Step 5).'
# Modify vcf so all chromosomes are prefixed with 'chr'
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' "$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.vcf > "$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.withChr.vcf

# At this point, the VCF, all.fasta and chain files are prefixed with 'chr' - so they are compatible.

mv "$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.withChr.vcf "$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.vcf

pigz --best -k "$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.vcf

echo '19. Converting VCF to be on maternal and paternal reference genome coordinate systems (Maternal step).'
# Remap VCF to parental genome coordinates using modified CrossMap (fixedBugs)
python ../software/CrossMap-0.2.3/bin/CrossMap.py vcf "$SAMPLE_ID"/"$SAMPLE_ID".maternal.edit.chain "$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.vcf "$SAMPLE_ID"/maternal/"$SAMPLE_ID".maternal.renamed.fa "$SAMPLE_ID"/maternal/"$SAMPLE_ID".maternal.vcf
echo '20. Converting VCF to be on maternal and paternal reference genome coordinate systems (Paternal step).'
python ../software/CrossMap-0.2.3/bin/CrossMap.py vcf "$SAMPLE_ID"/"$SAMPLE_ID".paternal.edit.chain "$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.vcf "$SAMPLE_ID"/paternal/"$SAMPLE_ID".paternal.renamed.fa "$SAMPLE_ID"/paternal/"$SAMPLE_ID".paternal.vcf

# This script swaps the REF and ALT alleles according to whether its the maternal or paternal haplotype.
echo '21. Correct Alleles in VCF to be parent specific.'
python ../software/PersonalGenomePipeline/haplotypeVCFAlleles.py "$SAMPLE_ID" "$SAMPLE_ID".maternal.vcf "$SAMPLE_ID".paternal.vcf
#
mv "$SAMPLE_ID"/maternal/"$SAMPLE_ID".maternal.vcf.maternal2.vcf "$SAMPLE_ID"/maternal/"$SAMPLE_ID".maternal.vcf
mv "$SAMPLE_ID"/paternal/"$SAMPLE_ID".paternal.vcf.paternal2.vcf "$SAMPLE_ID"/paternal/"$SAMPLE_ID".paternal.vcf

echo '22. running CreateSequenceDictionary.'
java -jar /home/centos/scratch/software/picard-tools-2.4.1/picard.jar CreateSequenceDictionary R=""$SAMPLE_ID"/paternal/$SAMPLE_ID".paternal.renamed.fa O="$SAMPLE_ID"/paternal/"$SAMPLE_ID".paternal.renamed.dict
java -jar /home/centos/scratch/software/picard-tools-2.4.1/picard.jar CreateSequenceDictionary R=""$SAMPLE_ID"/maternal/$SAMPLE_ID".maternal.renamed.fa O="$SAMPLE_ID"/maternal/"$SAMPLE_ID".maternal.renamed.dict
echo '23. running AddOrReplaceReadGroups.'
java -jar /home/centos/scratch/software/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="$SAMPLE_ID"/maternal/"$SAMPLE_ID".consensus.mat.filtered.sorted.bam O="$SAMPLE_ID"/maternal/"$SAMPLE_ID".consensus.mat.filtered.sorted.readGroup.bam  RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
java -jar /home/centos/scratch/software/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="$SAMPLE_ID"/paternal/"$SAMPLE_ID".consensus.pat.filtered.sorted.bam O="$SAMPLE_ID"/paternal/"$SAMPLE_ID".consensus.pat.filtered.sorted.readGroup.bam  RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20

# both say maternal for chromosome because maternal was used as a template BAM when selecting best reads.
echo '24. Harmonizing filtered BAMs'
~/scratch/software/samtools-1.3.1/samtools view -h "$SAMPLE_ID"/maternal/"$SAMPLE_ID".consensus.mat.filtered.sorted.readGroup.bam  | sed -e 's/_maternal//g' >> "$SAMPLE_ID"/maternal/"$SAMPLE_ID".consensus.mat.filtered.sorted2.bam
mv "$SAMPLE_ID"/maternal/"$SAMPLE_ID".consensus.mat.filtered.sorted2.bam "$SAMPLE_ID"/maternal/"$SAMPLE_ID".consensus.mat.filtered.sorted.readGroup.bam

~/scratch/software/samtools-1.3.1/samtools view -h "$SAMPLE_ID"/paternal/"$SAMPLE_ID".consensus.pat.filtered.sorted.readGroup.bam  | sed -e 's/_maternal//g' >> "$SAMPLE_ID"/paternal/"$SAMPLE_ID".consensus.pat.filtered.sorted2.bam
mv "$SAMPLE_ID"/paternal/"$SAMPLE_ID".consensus.pat.filtered.sorted2.bam "$SAMPLE_ID"/paternal/"$SAMPLE_ID".consensus.pat.filtered.sorted.readGroup.bam

echo '25. running ReorderSam.'
java -jar /home/centos/scratch/software/picard-tools-2.4.1/picard.jar ReorderSam I="$SAMPLE_ID"/paternal/"$SAMPLE_ID".consensus.pat.filtered.sorted.readGroup.bam O="$SAMPLE_ID"/paternal/"$SAMPLE_ID".consensus.pat.filtered.sorted.readGroup.bam2 R="$SAMPLE_ID"/paternal/"$SAMPLE_ID".paternal.renamed.fa ALLOW_CONTIG_LENGTH_DISCORDANCE=true
java -jar /home/centos/scratch/software/picard-tools-2.4.1/picard.jar ReorderSam I="$SAMPLE_ID"/maternal/"$SAMPLE_ID".consensus.mat.filtered.sorted.readGroup.bam O="$SAMPLE_ID"/maternal/"$SAMPLE_ID".consensus.mat.filtered.sorted.readGroup.bam2 R="$SAMPLE_ID"/maternal/"$SAMPLE_ID".maternal.renamed.fa ALLOW_CONTIG_LENGTH_DISCORDANCE=true

mv "$SAMPLE_ID"/paternal/"$SAMPLE_ID".consensus.pat.filtered.sorted.readGroup.bam2 "$SAMPLE_ID"/paternal/"$SAMPLE_ID".consensus.pat.filtered.sorted.readGroup.bam
mv "$SAMPLE_ID"/maternal/"$SAMPLE_ID".consensus.mat.filtered.sorted.readGroup.bam2 "$SAMPLE_ID"/maternal/"$SAMPLE_ID".consensus.mat.filtered.sorted.readGroup.bam

~/scratch/software/samtools-1.3.1/samtools index "$SAMPLE_ID"/maternal/"$SAMPLE_ID".consensus.mat.filtered.sorted.readGroup.bam
~/scratch/software/samtools-1.3.1/samtools index "$SAMPLE_ID"/paternal/"$SAMPLE_ID".consensus.pat.filtered.sorted.readGroup.bam

echo '26. Calculating ASE (Maternal).'
java -jar ~/scratch/software/GenomeAnalysisTK.jar -R "$SAMPLE_ID"/maternal/"$SAMPLE_ID".maternal.renamed.fa -T ASEReadCounter -o "$SAMPLE_ID"/maternal/"$SAMPLE_ID".ASE.mat.csv -I "$SAMPLE_ID"/maternal/"$SAMPLE_ID".consensus.mat.filtered.sorted.readGroup.bam -sites "$SAMPLE_ID"/maternal/"$SAMPLE_ID".maternal.vcf -dels -U ALLOW_N_CIGAR_READS -S SILENT
echo '27. Calculating ASE (Paternal).'
java -jar ~/scratch/software/GenomeAnalysisTK.jar -R "$SAMPLE_ID"/paternal/"$SAMPLE_ID".paternal.renamed.fa -T ASEReadCounter -o "$SAMPLE_ID"/paternal/"$SAMPLE_ID".ASE.pat.csv -I "$SAMPLE_ID"/paternal/"$SAMPLE_ID".consensus.pat.filtered.sorted.readGroup.bam -sites "$SAMPLE_ID"/paternal/"$SAMPLE_ID".paternal.vcf -dels -U ALLOW_N_CIGAR_READS -S SILENT

   #Currently, the ASE output is relative to the haplotype reference used (Maternal reference or paternal reference) 
   #i.e. REF in ASE.mat.csv = maternal allele etc.To assess REF Bias (which we should not expect due to using Personal genomes), 
   #and for easier interpretation, its better to have ASE in terms of universal REF / ALT. i.e. adding up all the REF/ALT counts 
   #for each haplotype, consolidated into a single ASE file per individual.

echo '28. Converting coordinates to reference genome.'
python ../software/PersonalGenomePipeline/ASERefCord.py "$SAMPLE_ID"/maternal/"$SAMPLE_ID".ASE.mat.csv "$SAMPLE_ID"/paternal/"$SAMPLE_ID".ASE.pat.csv "$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.vcf "$SAMPLE_ID"/maternal/"$SAMPLE_ID".maternal.vcf "$SAMPLE_ID"/paternal/"$SAMPLE_ID".paternal.vcf "$SAMPLE_ID".maternal.alleles.csv "$SAMPLE_ID".paternal.alleles.csv "$SAMPLE_ID"/maternal/"$SAMPLE_ID".maternal.ref.csv "$SAMPLE_ID"/paternal/"$SAMPLE_ID".paternal.ref.csv
#################################
echo '29. Combining haplotypes to produce file ASE file'
Rscript ../software/PersonalGenomePipeline/ASERefOut.R "$SAMPLE_ID"/maternal/"$SAMPLE_ID".maternal.ref.csv "$SAMPLE_ID"/paternal/"$SAMPLE_ID".paternal.ref.csv "$SAMPLE_ID"/"$SAMPLE_ID".ASE.csv

### END OF ASE PIPELINE ###

### START OF GENE-LEVEL COUNT GENERATION ###
echo '30. LiftingOver GTF (gencode.v19) (Maternal).'
../software/liftOver -gff gencode.v19.annotation.gtf "$SAMPLE_ID"/maternal.chain "$SAMPLE_ID"/"$SAMPLE_ID".mat.gtf "$SAMPLE_ID".not_lifted.txt -minMatch=0.00000001
echo '31. LiftingOver GTF (gencode.v19) (Paternal).'
../software/liftOver -gff gencode.v19.annotation.gtf "$SAMPLE_ID"/paternal.chain "$SAMPLE_ID"/"$SAMPLE_ID".pat.gtf "$SAMPLE_ID".not_lifted.txt -minMatch=0.00000001

# remove _parental string from the chromosome column in the liftedOver GTFs so it matches with the BAM files.
 sed -i "s/_maternal//g" "$SAMPLE_ID"/"$SAMPLE_ID".mat.gtf
 sed -i "s/_paternal//g" "$SAMPLE_ID"/"$SAMPLE_ID".pat.gtf

# Gene-level counts without multi-mapping - for each haplotype. 

# Haplotypic feature count (gene level counts)
echo '32. generating Gene-level counts (Maternal).'
../software/subread-1.5.0-p3-Linux-x86_64/bin/featureCounts -p -T 8 -a "$SAMPLE_ID"/"$SAMPLE_ID".mat.gtf -o "$SAMPLE_ID"/maternal/"$SAMPLE_ID".GeneCount_Mat.txt "$SAMPLE_ID"/maternal/"$SAMPLE_ID".consensus.mat.filtered.sorted.readGroup.bam
echo '33. generating Gene-level counts (Paternal).'
../software/subread-1.5.0-p3-Linux-x86_64/bin/featureCounts -p -T 8 -a "$SAMPLE_ID"/"$SAMPLE_ID".pat.gtf -o "$SAMPLE_ID"/paternal/"$SAMPLE_ID".GeneCount_Pat.txt "$SAMPLE_ID"/paternal/"$SAMPLE_ID".consensus.pat.filtered.sorted.readGroup.bam
echo '34. generating Gene-level counts (Reference).'
../software/subread-1.5.0-p3-Linux-x86_64/bin/featureCounts -p -T 8 -a gencode.v19.annotation.gtf -o "$SAMPLE_ID"/reference/"$SAMPLE_ID".GeneCount_Ref.txt "$SAMPLE_ID"/reference/"$SAMPLE_ID".filtered.bam

echo '35. Adding both haplotypes together to produce final gene-level count file'
Rscript ../software/PersonalGenomePipeline/AddHaploCounts.R "$SAMPLE_ID" "$SAMPLE_ID"/maternal/"$SAMPLE_ID".GeneCount_Mat.txt "$SAMPLE_ID"/paternal/"$SAMPLE_ID".GeneCount_Pat.txt "$SAMPLE_ID"/"$SAMPLE_ID".GeneCount.Final.txt
echo '36. Finished successfully.'


# clean up step to save space.
rm "$SAMPLE_ID"/maternal/chr*_"$SAMPLE_ID"_maternal.fa
rm "$SAMPLE_ID"/paternal/chr*_"$SAMPLE_ID"_paternal.fa

mkdir "$SAMPLE_ID"/map_files
mv "$SAMPLE_ID"/*.map map_files/

rm "$SAMPLE_ID"/maternal/SA
rm "$SAMPLE_ID"/maternal/SAindex
rm "$SAMPLE_ID"/maternal/Genome

rm "$SAMPLE_ID"/paternal/SA
rm "$SAMPLE_ID"/paternal/SAindex
rm "$SAMPLE_ID"/paternal/Genome

rm "$SAMPLE_ID"/"$SAMPLE_ID"_1.fastq
rm "$SAMPLE_ID"/"$SAMPLE_ID"_2.fastq

" > $SAMPLE_ID.both.sh
fi

done < ${ALL_SAMPLES}
