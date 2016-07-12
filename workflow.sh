# USAGE:   sh SCRIPT.sh ALL_SAMPLEs.txt UK10K_SAMPLES.txt
# This should all be excuted from the /pvfs2/ directory.
THREAD_NO=$3
ALL_SAMPLES=$1
UK10K_SAMPLES=$2
FML='"chr"'
VAR="'{if(\$0 !~ /^#/) print $FML\$0; else print \$0}'"

while read line ; do

## check if line exist in CHECK_FILE; then assign result to variable
SAMPLE_ID=$(grep "^${line}$" ${UK10K_SAMPLES})


## if variable is blank (meaning ALL_SAMPLES line not found in UK10K)
## print 'false' and exit
if [[ -z $SAMPLE_ID ]] ; then

echo '#!/bin/bash 
# 
#SBATCH -N 1 
#SBATCH --mail-type=END
#SBATCH --mail-user=craig.glastonbury@kcl.ac.uk
# number of nodes 
#SBATCH -n '$THREAD_NO' 

echo "Step 1. Sorting BAM file"
/media/shared_data/software/samtools-1.3.1/samtools sort -@ 6 -m 10G -n '$line'/'$line'_sorted.bam -o '$line'/'$line'_sorted.bam.sorted
echo "Step 2. BAM >> fastq"
/media/shared_data/software/bedtools2/bin/bedtools bamtofastq -i '$line'/'$line'_sorted.bam.sorted -fq '$line'/'$line'.f2.fq -fq2 '$line'/'$line'.f1.fq
echo "Step 3. Trimming Adaptor seqs"
/media/shared_data/software/trim_galore_zip/trim_galore -stringency 5 -q 1 -o '$line' --paired '$line'/'$line'.f1.fq '$line'/'$line'.f2.fq
echo "Step 4. Trimming polyA+ seqs"
perl /media/shared_data/software/prinseq-lite-0.20.4/prinseq-lite.pl -fastq '$line'/'$line'.f1_val_1.fq -fastq2 '$line'/'$line'.f2_val_2.fq -out_good '$line'/'$line' -trim_tail_left 5 -trim_tail_right 5 -min_len 20
rm '$line'/'$line'.f2_val_2.fq
rm '$line'/'$line'.f1_val_1.fq
rm '$line'/'$line'.f2.fq_trimming_report.txt
rm '$line'/'$line'.f1.fq_trimming_report.txt
rm '$line'/'$line'_1_singletons.fastq
rm '$line'/'$line'_2_singletons.fastq
rm '$line'/'$line'_prinseq_bad*
rm '$line'/'$line'.f1.fq
rm '$line'/'$line'.f2.fq
rm '$line'/'$line'_sorted.bam.sorted
echo "Step 5. Genome Alignment (STAR)"
/media/shared_data/software/STAR/bin/Linux_x86_64/STAR --runThreadN '$THREAD_NO' --runMode alignReads --readFilesIn '$line'/'$line'_1.fastq '$line'/'$line'_2.fastq --genomeDir hg19 --outSAMstrandField intronMotif --outFilterMultimapNmax 30 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 15 --outMultimapperOrder Random --outSAMunmapped Within --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --sjdbOverhang 48 --outFilterMismatchNmax 6 --outSAMattributes NH nM NM MD HI --outSAMattrRGline  ID:'$line'_maternal PU:Illumina PL:Illumina LB:'$line'_maternal SM:'$line'_maternal CN:Seq_centre --outSAMtype BAM SortedByCoordinate --outFileNamePrefix '$line'_ref.

mv '$line'_ref.Aligned.sortedByCoord.out.bam '$line'/reference/'$line'_ref.Aligned.sortedByCoord.out.bam
mv '$line'_ref.Chimeric.out.junction '$line'/reference/'$line'_ref.Chimeric.out.junction
mv '$line'_ref.Chimeric.out.sam '$line'/reference/'$line'_ref.Chimeric.out.sam
mv '$line'_ref.Log.final.out '$line'/reference/'$line'_ref.Log.final.out
mv '$line'_ref.SJ.out.tab '$line'/reference/'$line'_ref.SJ.out.tab
rm '$line'_ref.Log.out
rm '$line'_ref.Log.progress.out

echo "Step 6. Filtering aligned BAM"
/media/shared_data/software/samtools-1.3.1/samtools view -@ 6 -b -F4 -q 30 '$line'/reference/'$line'_ref.Aligned.sortedByCoord.out.bam -o '$line'/reference/'$line'.filtered.bam

echo "Step 7. Calculating gene counts"
/media/shared_data/software/subread-1.5.0-p3-Linux-x86_64/bin/featureCounts -p -T 8 -a gencode.v19.annotation.gtf -o '$line'/reference/'$line'.GeneCount_Ref.txt '$line'/reference/'$line'.filtered.bam

rm '$line'/'$line'_sorted.bam.sorted
rm '$line'/'$line'_sorted.bam 
rm -r '$line'_ref._STARtmp
rm '$line'/'$line'_1.fastq
rm '$line'/'$line'_2.fastq

rm '$line'/reference/'$line'.filtered.bam' > $line.refOnly.sh

else
echo '#!/bin/bash 
# 
#SBATCH -N 1 
# number of nodes 
#SBATCH -n '$THREAD_NO' 
#SBATCH --mail-type=END
#SBATCH --mail-user=craig.glastonbury@kcl.ac.uk

# Only edit the genome using variants with a GQ score > 30.

echo "Step 1. Filtering VCF for GQ < 30"
python /media/shared_data/software/PersonalGenomePipeline/filterVariantGQ.py '$SAMPLE_ID'

mv '$SAMPLE_ID'/'$SAMPLE_ID'.vcf.gz2 '$SAMPLE_ID'/'$SAMPLE_ID'.vcf

rm '$SAMPLE_ID'/'$SAMPLE_ID'.vcf.gz

pigz --best -k '$SAMPLE_ID'/'$SAMPLE_ID'.vcf

echo "Step 2. Personalising genomes with vcf2diploid"
java -jar /media/shared_data/software/vcf2diploid_v0.2.6a/vcf2diploid.jar -id '$SAMPLE_ID' -chr hg19/hg19.fa -vcf '$SAMPLE_ID'/'$SAMPLE_ID'.vcf.gz -outDir '$SAMPLE_ID'

mv '$SAMPLE_ID'/*_'$SAMPLE_ID'_maternal.fa '$SAMPLE_ID'/maternal/
mv '$SAMPLE_ID'/*_'$SAMPLE_ID'_paternal.fa '$SAMPLE_ID'/paternal/

# Make BAMS sorted by read name

echo "Step 3. Sorting BAM by read name"
/media/shared_data/software/samtools-1.3.1/samtools sort -@ 6 -m 10G -n '$SAMPLE_ID'/'$SAMPLE_ID'_sorted.bam -o '$SAMPLE_ID'/'$SAMPLE_ID'_sorted.bam.sorted

# convert BAMS to fastq 
echo "Step 4. BAM >> fastq"
/media/shared_data/software/bedtools2/bin/bedtools bamtofastq -i '$SAMPLE_ID'/'$SAMPLE_ID'_sorted.bam.sorted -fq '$SAMPLE_ID'/'$SAMPLE_ID'.f2.fq -fq2 '$SAMPLE_ID'/'$SAMPLE_ID'.f1.fq

# Add check to see whether both ref and personal alignments should be done, or just one (i.e. individuals not in UK10K)
echo "Step 5. Trimming adaptor seqs"
/media/shared_data/software/trim_galore_zip/trim_galore -stringency 5 -q 1 -o '$SAMPLE_ID' --phred33 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --paired '$SAMPLE_ID'/'$SAMPLE_ID'.f1.fq '$SAMPLE_ID'/'$SAMPLE_ID'.f2.fq

echo "Step 6. Trimming PolyA+ seqs"
perl /media/shared_data/software/prinseq-lite-0.20.4/prinseq-lite.pl -fastq '$SAMPLE_ID'/'$SAMPLE_ID'.f1_val_1.fq -fastq2 '$SAMPLE_ID'/'$SAMPLE_ID'.f2_val_2.fq -out_good '$SAMPLE_ID'/'$SAMPLE_ID' -trim_tail_left 5 -trim_tail_right 5 -min_len 20

rm '$SAMPLE_ID'/'$SAMPLE_ID'.f2_val_2.fq
rm '$SAMPLE_ID'/'$SAMPLE_ID'.f1_val_1.fq
rm '$SAMPLE_ID'/'$SAMPLE_ID'.f2.fq_trimming_report.txt
rm '$SAMPLE_ID'/'$SAMPLE_ID'.f1.fq_trimming_report.txt
rm '$SAMPLE_ID'/'$SAMPLE_ID'_1_singletons.fastq
rm '$SAMPLE_ID'/'$SAMPLE_ID'_2_singletons.fastq
rm '$SAMPLE_ID'/'$SAMPLE_ID'_prinseq_bad*
rm '$SAMPLE_ID'/'$SAMPLE_ID'.f1.fq
rm '$SAMPLE_ID'/'$SAMPLE_ID'.f2.fq
rm '$SAMPLE_ID'/'$SAMPLE_ID'_sorted.bam.sorted # Ensure this is the name after replacing BAM_IDs with TwinUK IDs.
# Maternal genome generation (suffix arrays etc)

echo "Step 7. Generating Maternal genome Suffix Array (STAR)"
/media/shared_data/software/STAR/bin/Linux_x86_64/STAR --runThreadN '$THREAD_NO' --runMode genomeGenerate --genomeDir '$SAMPLE_ID'/maternal --genomeFastaFiles '$SAMPLE_ID'/maternal/chr*_'$SAMPLE_ID'_maternal.fa

#Paternal genome generation (suffix arrays etc)
echo "Step 8. Generating Paternal genome Suffix Array (STAR)"
/media/shared_data/software/STAR/bin/Linux_x86_64/STAR --runThreadN '$THREAD_NO' --runMode genomeGenerate --genomeDir '$SAMPLE_ID'/paternal --genomeFastaFiles '$SAMPLE_ID'/paternal/chr*_'$SAMPLE_ID'_paternal.fa

#Align paternal
echo "Step 9. Aligning to Paternal genome (STAR)"
/media/shared_data/software/STAR/bin/Linux_x86_64/STAR --runThreadN '$THREAD_NO' --runMode alignReads --readFilesIn '$SAMPLE_ID'/'$SAMPLE_ID'_1.fastq '$SAMPLE_ID'/'$SAMPLE_ID'_2.fastq --genomeDir '$SAMPLE_ID'/paternal --outSAMstrandField intronMotif --outFilterMultimapNmax 30 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 15 --outMultimapperOrder Random --outSAMunmapped Within --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --sjdbOverhang 48 --outFilterMismatchNmax 6 --outSAMattributes NH nM NM MD HI --outSAMattrRGline  ID:'$SAMPLE_ID'_maternal PU:Illumina PL:Illumina LB:'$SAMPLE_ID'_maternal SM:'$SAMPLE_ID'_maternal CN:Seq_centre --outSAMtype BAM SortedByCoordinate --outFileNamePrefix '$SAMPLE_ID'_pat.
mv '$SAMPLE_ID'_pat.Aligned.sortedByCoord.out.bam '$SAMPLE_ID'/paternal/'$SAMPLE_ID'_pat.Aligned.sortedByCoord.out.bam
mv '$SAMPLE_ID'_pat.Chimeric.out.junction '$SAMPLE_ID'/paternal/'$SAMPLE_ID'_pat.Chimeric.out.junction
mv '$SAMPLE_ID'_pat.Chimeric.out.sam '$SAMPLE_ID'/paternal/'$SAMPLE_ID'_pat.Chimeric.out.sam
mv '$SAMPLE_ID'_pat.Log.final.out '$SAMPLE_ID'/paternal/'$SAMPLE_ID'_pat.Log.final.out
mv '$SAMPLE_ID'_pat.SJ.out.tab '$SAMPLE_ID'/paternal/'$SAMPLE_ID'_pat.SJ.out.tab
rm '$SAMPLE_ID'_pat.Log.out
rm '$SAMPLE_ID'_pat.Log.progress.out

echo "Step 10. Aligning to Maternal genome (STAR)"

/media/shared_data/software/STAR/bin/Linux_x86_64/STAR --runThreadN '$THREAD_NO' --runMode alignReads --readFilesIn '$SAMPLE_ID'/'$SAMPLE_ID'_1.fastq '$SAMPLE_ID'/'$SAMPLE_ID'_2.fastq --genomeDir '$SAMPLE_ID'/maternal --outSAMstrandField intronMotif --outFilterMultimapNmax 30 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 15 --outMultimapperOrder Random --outSAMunmapped Within --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --sjdbOverhang 48 --outFilterMismatchNmax 6 --outSAMattributes NH nM NM MD HI --outSAMattrRGline  ID:'$SAMPLE_ID'_maternal PU:Illumina PL:Illumina LB:'$SAMPLE_ID'_maternal SM:'$SAMPLE_ID'_maternal CN:Seq_centre --outSAMtype BAM SortedByCoordinate --outFileNamePrefix '$SAMPLE_ID'_mat.
mv '$SAMPLE_ID'_mat.Aligned.sortedByCoord.out.bam '$SAMPLE_ID'/maternal/'$SAMPLE_ID'_mat.Aligned.sortedByCoord.out.bam
mv '$SAMPLE_ID'_mat.Chimeric.out.junction '$SAMPLE_ID'/maternal/'$SAMPLE_ID'_mat.Chimeric.out.junction
mv '$SAMPLE_ID'_mat.Chimeric.out.sam '$SAMPLE_ID'/maternal/'$SAMPLE_ID'_mat.Chimeric.out.sam
mv '$SAMPLE_ID'_mat.Log.final.out '$SAMPLE_ID'/maternal/'$SAMPLE_ID'_mat.Log.final.out
mv '$SAMPLE_ID'_mat.SJ.out.tab '$SAMPLE_ID'/maternal/'$SAMPLE_ID'_mat.SJ.out.tab
rm '$SAMPLE_ID'_mat.Log.out
rm '$SAMPLE_ID'_mat.Log.progress.out

# Add standard reference alignment too.

echo "Step 11. Aligning to Universal Reference genome (STAR)"

/media/shared_data/software/STAR/bin/Linux_x86_64/STAR --runThreadN '$THREAD_NO' --runMode alignReads --readFilesIn '$SAMPLE_ID'/'$SAMPLE_ID'_1.fastq '$SAMPLE_ID'/'$SAMPLE_ID'_2.fastq --genomeDir hg19 --outSAMstrandField intronMotif --outFilterMultimapNmax 30 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 15 --outMultimapperOrder Random --outSAMunmapped Within --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --sjdbOverhang 48 --outFilterMismatchNmax 6 --outSAMattributes NH nM NM MD HI --outSAMattrRGline  ID:'$SAMPLE_ID'_ref PU:Illumina PL:Illumina LB:'$SAMPLE_ID'_ref SM:'$SAMPLE_ID'_ref CN:Seq_centre --outSAMtype BAM SortedByCoordinate --outFileNamePrefix '$SAMPLE_ID'_ref.

mv '$SAMPLE_ID'_ref.Aligned.sortedByCoord.out.bam '$SAMPLE_ID'/reference/'$SAMPLE_ID'_ref.Aligned.sortedByCoord.out.bam
mv '$SAMPLE_ID'_ref.Chimeric.out.junction '$SAMPLE_ID'/reference/'$SAMPLE_ID'_ref.Chimeric.out.junction
mv '$SAMPLE_ID'_ref.Chimeric.out.sam '$SAMPLE_ID'/reference/'$SAMPLE_ID'_ref.Chimeric.out.sam
mv '$SAMPLE_ID'_ref.Log.final.out '$SAMPLE_ID'/reference/'$SAMPLE_ID'_ref.Log.final.out
mv '$SAMPLE_ID'_ref.SJ.out.tab '$SAMPLE_ID'/reference/'$SAMPLE_ID'_ref.SJ.out.tab
rm '$SAMPLE_ID'_ref.Log.out
rm '$SAMPLE_ID'_ref.Log.progress.out
########################################
#
rm '$SAMPLE_ID'/maternal/chrLength.txt
rm '$SAMPLE_ID'/maternal/chrNameLength.txt
rm '$SAMPLE_ID'/maternal/chrStart.txt
rm '$SAMPLE_ID'/maternal/Genome 
rm '$SAMPLE_ID'/maternal/genomeParameters.txt
rm '$SAMPLE_ID'/maternal/SA
rm '$SAMPLE_ID'/maternal/SAindex                  

rm '$SAMPLE_ID'/paternal/chrLength.txt
rm '$SAMPLE_ID'/paternal/chrNameLength.txt
rm '$SAMPLE_ID'/paternal/chrStart.txt
rm '$SAMPLE_ID'/paternal/Genome
rm '$SAMPLE_ID'/paternal/genomeParameters.txt
rm '$SAMPLE_ID'/paternal/SA                  
rm '$SAMPLE_ID'/paternal/SAindex


echo "Step 12. Filtering BAM file"
/media/shared_data/software/samtools-1.3.1/samtools view -@ 6 -b -F4 -q 30 '$SAMPLE_ID'/reference/'$SAMPLE_ID'_ref.Aligned.sortedByCoord.out.bam -o '$SAMPLE_ID'/reference/'$SAMPLE_ID'.filtered.bam

## Produce consensus bams, in which the best read per haplotype is selected.

echo "Step 13. Selecting best read alignment per haplotype (seperateBAMS.py)"
python /media/shared_data/software/PersonalGenomePipeline/pipeline.seperateBAMs.py '$SAMPLE_ID'
#
#for all variants without an rsid - assign them chrpos.
echo "Step 14. Assigning variants with no rsid to: chr-pos"
python /media/shared_data/software/PersonalGenomePipeline/nameVariants.py '$SAMPLE_ID' '$SAMPLE_ID'.hets.phased.vcf.gz '$SAMPLE_ID'.hets.phased.vcf2.gz

mv '$SAMPLE_ID'/'$SAMPLE_ID'.hets.phased.vcf2.gz '$SAMPLE_ID'/'$SAMPLE_ID'.hets.phased.vcf

# Sort the VCFs else picard/GATK throws a fit.
echo "Step 15. Sorting VCF by chromosome"
/media/shared_data/software/vcftools/src/perl/vcf-sort -c '$SAMPLE_ID'/'$SAMPLE_ID'.hets.phased.vcf > '$SAMPLE_ID'/'$SAMPLE_ID'.hets.GATK.sorted.vcf

#########################################
# run BASH script to concatenate reference genomes, and rename chain 1_maternal > 1

echo "Step 16. Consolidating nomenclature across files (so all files are prefixed with chr)"
cp '$SAMPLE_ID'/paternal.chain '$SAMPLE_ID'/'$SAMPLE_ID'.paternal.edit.chain
cp '$SAMPLE_ID'/maternal.chain '$SAMPLE_ID'/'$SAMPLE_ID'.maternal.edit.chain
for i in {1..22}
do
	sed -i "s/"$i"_maternal/"$i"/g" '$SAMPLE_ID'/'$SAMPLE_ID'.maternal.edit.chain
	sed -i "s/"$i"_paternal/"$i"/g" '$SAMPLE_ID'/'$SAMPLE_ID'.paternal.edit.chain
	cat '$SAMPLE_ID'/paternal/chr"$i"_'$SAMPLE_ID'_paternal.fa >> '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.all_paternal.fa
	cat '$SAMPLE_ID'/maternal/chr"$i"_'$SAMPLE_ID'_maternal.fa >> '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.all_maternal.fa

done

sed -i 's/X_paternal/X/g' '$SAMPLE_ID'/'$SAMPLE_ID'.paternal.edit.chain
#sed -i 's/Y_paternal/Y/g' '$SAMPLE_ID'/'$SAMPLE_ID'.paternal.edit.chain
sed -i 's/M_paternal/MT/g' '$SAMPLE_ID'/'$SAMPLE_ID'.paternal.edit.chain

sed -i 's/X_maternal/X/g' '$SAMPLE_ID'/'$SAMPLE_ID'.maternal.edit.chain
#sed -i 's/Y_maternal/Y/g' '$SAMPLE_ID'/'$SAMPLE_ID'.maternal.edit.chain
sed -i 's/M_maternal/MT/g' '$SAMPLE_ID'/'$SAMPLE_ID'.maternal.edit.chain

cat '$SAMPLE_ID'/paternal/chrX_'$SAMPLE_ID'_paternal.fa >> '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.all_paternal.fa
#cat '$SAMPLE_ID'/paternal/chrY_'$SAMPLE_ID'_paternal.fa >> '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.all_paternal.fa
cat '$SAMPLE_ID'/paternal/chrM_'$SAMPLE_ID'_paternal.fa >> '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.all_paternal.fa

cat '$SAMPLE_ID'/maternal/chrX_'$SAMPLE_ID'_maternal.fa >> '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.all_maternal.fa
#cat '$SAMPLE_ID'/maternal/chrY_'$SAMPLE_ID'_maternal.fa >> '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.all_maternal.fa
cat '$SAMPLE_ID'/maternal/chrM_'$SAMPLE_ID'_maternal.fa >> '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.all_maternal.fa

#################################

# rename maternal and paternal fasta files with correct scaffold/chr name.

python /media/shared_data/software/PersonalGenomePipeline/renameFaChr.py '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.all_maternal.fa '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.maternal.renamed.fa '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.all_paternal.fa '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.paternal.renamed.fa

#pigz -d '$SAMPLE_ID'/'$SAMPLE_ID'.hets.GATK.sorted.vcf.gz

# Modify vcf so all chromosomes are prefixed with 'chr'
awk '$VAR' '$SAMPLE_ID'/'$SAMPLE_ID'.hets.GATK.sorted.vcf > '$SAMPLE_ID'/'$SAMPLE_ID'.hets.GATK.sorted.withChr.vcf

# At this point, the VCF, all.fasta and chain files are prefixed with 'chr' - so they are compatible.

mv '$SAMPLE_ID'/'$SAMPLE_ID'.hets.GATK.sorted.withChr.vcf '$SAMPLE_ID'/'$SAMPLE_ID'.hets.GATK.sorted.vcf

# Remap VCF to parental genome coordinates using modified CrossMap (fixedBugs)

echo "Step 17. Converting VCF to Maternal coordinate system for downstream ASE"
python /media/shared_data/software/CrossMap-0.2.3/bin/CrossMap.py vcf '$SAMPLE_ID'/'$SAMPLE_ID'.maternal.edit.chain '$SAMPLE_ID'/'$SAMPLE_ID'.hets.GATK.sorted.vcf '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.maternal.renamed.fa '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.maternal.vcf

echo "Step 18. Converting VCF to Paternal coordinate system for downstream ASE"
python /media/shared_data/software/CrossMap-0.2.3/bin/CrossMap.py vcf '$SAMPLE_ID'/'$SAMPLE_ID'.paternal.edit.chain '$SAMPLE_ID'/'$SAMPLE_ID'.hets.GATK.sorted.vcf '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.paternal.renamed.fa '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.paternal.vcf

# This script swaps the REF and ALT alleles according to whether its the maternal or paternal haplotype.
echo "Step 19. Swapping Alleles to reflect reference according to haplotype"
python /media/shared_data/software/PersonalGenomePipeline/haplotypeVCFAlleles.py '$SAMPLE_ID' '$SAMPLE_ID'.maternal.vcf '$SAMPLE_ID'.paternal.vcf
#
mv '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.maternal.vcf.maternal2.vcf '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.maternal.vcf
mv '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.paternal.vcf.paternal2.vcf '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.paternal.vcf

echo "Step 20. Creating SequenceDict (picard)"
java -jar /media/shared_data/software/picard-tools-2.4.1/picard.jar CreateSequenceDictionary R='$SAMPLE_ID'/paternal/'$SAMPLE_ID'.paternal.renamed.fa O='$SAMPLE_ID'/paternal/'$SAMPLE_ID'.paternal.renamed.dict
java -jar /media/shared_data/software/picard-tools-2.4.1/picard.jar CreateSequenceDictionary R='$SAMPLE_ID'/maternal/'$SAMPLE_ID'.maternal.renamed.fa O='$SAMPLE_ID'/maternal/'$SAMPLE_ID'.maternal.renamed.dict

echo "Step 21. Adding ReadGroups (picard)"
java -jar /media/shared_data/software/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I='$SAMPLE_ID'/maternal/'$SAMPLE_ID'.consensus.mat.filtered.sorted.bam O='$SAMPLE_ID'/maternal/'$SAMPLE_ID'.consensus.mat.filtered.sorted.readGroup.bam  RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
java -jar /media/shared_data/software/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I='$SAMPLE_ID'/paternal/'$SAMPLE_ID'.consensus.pat.filtered.sorted.bam O='$SAMPLE_ID'/paternal/'$SAMPLE_ID'.consensus.pat.filtered.sorted.readGroup.bam  RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20

# both say maternal for chromosome because maternal was used as a template BAM when selecting best reads.
/media/shared_data/software/samtools-1.3.1/samtools view -@ 6 -h '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.consensus.mat.filtered.sorted.readGroup.bam  | sed -e 's/_maternal//g' >> '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.consensus.mat.filtered.sorted2.bam
mv '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.consensus.mat.filtered.sorted2.bam '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.consensus.mat.filtered.sorted.readGroup.bam

/media/shared_data/software/samtools-1.3.1/samtools view -@ 6 -h '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.consensus.pat.filtered.sorted.readGroup.bam  | sed -e 's/_maternal//g' >> '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.consensus.pat.filtered.sorted2.bam
mv '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.consensus.pat.filtered.sorted2.bam '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.consensus.pat.filtered.sorted.readGroup.bam

echo "Step 22. Ordering SAM files"
java -jar /media/shared_data/software/picard-tools-2.4.1/picard.jar ReorderSam I='$SAMPLE_ID'/paternal/'$SAMPLE_ID'.consensus.pat.filtered.sorted.readGroup.bam O='$SAMPLE_ID'/paternal/'$SAMPLE_ID'.consensus.pat.filtered.sorted.readGroup.bam2 R='$SAMPLE_ID'/paternal/'$SAMPLE_ID'.paternal.renamed.fa ALLOW_CONTIG_LENGTH_DISCORDANCE=true
java -jar /media/shared_data/software/picard-tools-2.4.1/picard.jar ReorderSam I='$SAMPLE_ID'/maternal/'$SAMPLE_ID'.consensus.mat.filtered.sorted.readGroup.bam O='$SAMPLE_ID'/maternal/'$SAMPLE_ID'.consensus.mat.filtered.sorted.readGroup.bam2 R='$SAMPLE_ID'/maternal/'$SAMPLE_ID'.maternal.renamed.fa ALLOW_CONTIG_LENGTH_DISCORDANCE=true

mv '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.consensus.pat.filtered.sorted.readGroup.bam2 '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.consensus.pat.filtered.sorted.readGroup.bam
mv '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.consensus.mat.filtered.sorted.readGroup.bam2 '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.consensus.mat.filtered.sorted.readGroup.bam

echo "Step 23. Indexing SAM files"
/media/shared_data/software/samtools-1.3.1/samtools index '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.consensus.mat.filtered.sorted.readGroup.bam
/media/shared_data/software/samtools-1.3.1/samtools index '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.consensus.pat.filtered.sorted.readGroup.bam

echo "Step 24. Performing ASEReadCounter - Maternal"
java -jar /media/shared_data/software/GenomeAnalysisTK.jar -R '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.maternal.renamed.fa -T ASEReadCounter -o '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.ASE.mat.csv -I '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.consensus.mat.filtered.sorted.readGroup.bam -sites '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.maternal.vcf -dels -U ALLOW_N_CIGAR_READS -S SILENT

echo "Step 25. Performing ASEReadCounter - Paternal"
java -jar /media/shared_data/software/GenomeAnalysisTK.jar -R '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.paternal.renamed.fa -T ASEReadCounter -o '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.ASE.pat.csv -I '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.consensus.pat.filtered.sorted.readGroup.bam -sites '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.paternal.vcf -dels -U ALLOW_N_CIGAR_READS -S SILENT

   #Currently, the ASE output is relative to the haplotype reference used (Maternal reference or paternal reference) 
   #i.e. REF in ASE.mat.csv = maternal allele etc.To assess REF Bias (which we should not expect due to using Personal genomes), 
   #and for easier interpretation, its better to have ASE in terms of universal REF / ALT. i.e. adding up all the REF/ALT counts 
   #for each haplotype, consolidated into a single ASE file per individual.

echo "Step 26. Converting Maternal and Paternal coordinates onto reference-base"
python /media/shared_data/software/PersonalGenomePipeline/ASERefCord.py '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.ASE.mat.csv '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.ASE.pat.csv '$SAMPLE_ID'/'$SAMPLE_ID'.hets.GATK.sorted.vcf '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.maternal.vcf '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.paternal.vcf '$SAMPLE_ID'.maternal.alleles.csv '$SAMPLE_ID'.paternal.alleles.csv '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.maternal.ref.csv '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.paternal.ref.csv '$SAMPLE_ID'
#################################

echo "Step 27. Sorting BAMs by position"
/media/shared_data/software/samtools-1.3.1/samtools sort -@ 6 -m 10G 1152/paternal/'$SAMPLE_ID'_pat.Aligned.sortedByCoord.out.bam -o '$SAMPLE_ID'/paternal/'$SAMPLE_ID'_pat.Aligned.sortedByCoord.out.bam.sorted
/media/shared_data/software/samtools-1.3.1/samtools sort -@ 6 -m 10G 1152/maternal/'$SAMPLE_ID'_mat.Aligned.sortedByCoord.out.bam -o '$SAMPLE_ID'/maternal/'$SAMPLE_ID'_mat.Aligned.sortedByCoord.out.bam.sorted

mv '$SAMPLE_ID'/maternal/'$SAMPLE_ID'_mat.Aligned.sortedByCoord.out.bam.sorted '$SAMPLE_ID'/maternal/'$SAMPLE_ID'_mat.Aligned.sortedByCoord.out.bam
mv '$SAMPLE_ID'/paternal/'$SAMPLE_ID'_pat.Aligned.sortedByCoord.out.bam.sorted '$SAMPLE_ID'/paternal/'$SAMPLE_ID'_pat.Aligned.sortedByCoord.out.bam

echo "Step 28. Indexing BAMs "
/media/shared_data/software/samtools-1.3.1/samtools index '$SAMPLE_ID'/maternal/'$SAMPLE_ID'_mat.Aligned.sortedByCoord.out.bam
/media/shared_data/software/samtools-1.3.1/samtools index '$SAMPLE_ID'/paternal/'$SAMPLE_ID'_pat.Aligned.sortedByCoord.out.bam

# Add total number of reads covering heterozgyous SNPs + degree of multi-mapping reads.
echo "Step 29. Quantifying total read overlapping het sites + degree of haplotype specific multi-mapping"
python /media/shared_data/software/PersonalGenomePipeline/quantMultiMapping.py '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.maternal.ref.csv '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.paternal.ref.csv '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.ASE.mat.ref.multi.txt '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.ASE.pat.ref.multi.txt '$SAMPLE_ID'

echo "Step 30. Producing file ASE output"
Rscript /media/shared_data/software/PersonalGenomePipeline/ASERefOut.R '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.ASE.mat.ref.multi.txt '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.ASE.pat.ref.multi.txt '$SAMPLE_ID'/'$SAMPLE_ID'.ASE.csv

### END OF ASE PIPELINE ###

### START OF GENE-LEVEL COUNT GENERATION ###
echo "Step 31. Lifting over GENCODEv19 GTF to Maternal coordinates"
/media/shared_data/software/liftOver -gff gencode.v19.annotation.gtf '$SAMPLE_ID'/maternal.chain '$SAMPLE_ID'/'$SAMPLE_ID'.mat.gtf '$SAMPLE_ID'.not_lifted.txt -minMatch=0.00000001
echo "Step 32. Lifting over GENCODEv19 GTF to Paternal coordinates"
/media/shared_data/software/liftOver -gff gencode.v19.annotation.gtf '$SAMPLE_ID'/paternal.chain '$SAMPLE_ID'/'$SAMPLE_ID'.pat.gtf '$SAMPLE_ID'.not_lifted.txt -minMatch=0.00000001

# remove _parental string from the chromosome column in the liftedOver GTFs so it matches with the BAM files.
 sed -i "s/_maternal//g" '$SAMPLE_ID'/'$SAMPLE_ID'.mat.gtf
 sed -i "s/_paternal//g" '$SAMPLE_ID'/'$SAMPLE_ID'.pat.gtf

# Gene-level counts without multi-mapping - for each haplotype. 

# Haplotypic feature count (gene level counts)
echo "Step 33. Calculating Gene level counts - Maternal"
/media/shared_data/software/subread-1.5.0-p3-Linux-x86_64/bin/featureCounts -p -T 8 -a '$SAMPLE_ID'/'$SAMPLE_ID'.mat.gtf -o '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.GeneCount_Mat.txt '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.consensus.mat.filtered.sorted.readGroup.bam
echo "Step 34. Calculating Gene level counts - Paternal"
/media/shared_data/software/subread-1.5.0-p3-Linux-x86_64/bin/featureCounts -p -T 8 -a '$SAMPLE_ID'/'$SAMPLE_ID'.pat.gtf -o '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.GeneCount_Pat.txt '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.consensus.pat.filtered.sorted.readGroup.bam
echo "Step 35. Calculating Gene level counts - Universal Reference"
/media/shared_data/software/subread-1.5.0-p3-Linux-x86_64/bin/featureCounts -p -T 8 -a gencode.v19.annotation.gtf -o '$SAMPLE_ID'/reference/'$SAMPLE_ID'.GeneCount_Ref.txt '$SAMPLE_ID'/reference/'$SAMPLE_ID'.filtered.bam

echo "Step 36. Combining Maternal and Paternal gene-level counts"
Rscript /media/shared_data/software/PersonalGenomePipeline/AddHaploCounts.R '$SAMPLE_ID' '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.GeneCount_Mat.txt '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.GeneCount_Pat.txt '$SAMPLE_ID'/'$SAMPLE_ID'.GeneCount.Final.txt

rm -r '$SAMPLE_ID'_pat._STARtmp
rm -r '$SAMPLE_ID'_ref._STARtmp
rm -r '$SAMPLE_ID'_mat._STARtmp
rm '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.consensus.mat.filtered.sorted.readGroup.bam
rm '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.consensus.pat.filtered.sorted.readGroup.bam
rm '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.paternal.renamed.fa
rm '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.maternal.renamed.fa


rm '$SAMPLE_ID'/'$SAMPLE_ID'.pat.gtf
rm '$SAMPLE_ID'/'$SAMPLE_ID'.mat.gtf

echo "Step 37. Successfully completed."

rm '$SAMPLE_ID'.not_lifted.txt ' > $line.both.sh
fi
done < ${ALL_SAMPLES}
