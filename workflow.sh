# This should all be excuted from the /data/ directory.
THREAD_NO=8
TISSUE='adipose'

java -jar ../software/vcf2diploid_v0.2.6a/vcf2diploid.jar -id "$SAMPLE_ID" -chr hg19/hg19.fa -vcf "$SAMPLE_ID"/"$SAMPLE_ID".vcf.gz -outDir $SAMPLE_ID

mv "$TISSUE"/"$SAMPLE_ID"/*_"$SAMPLE_ID"_maternal.fa "$TISSUE"/"$SAMPLE_ID"/maternal/
mv "$TISSUE"/"$SAMPLE_ID"/*_"$SAMPLE_ID"_paternal.fa "$TISSUE"/"$SAMPLE_ID"/paternal/

# Make BAMS sorted by read name
../software/samtools-1.3.1/samtools sort -n "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID"_sorted.bam -o "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID"_sorted.bam.sorted

# convert BAMS to fastq 
../software/bedtools2/bin/bedtools bamtofastq -i "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID"_sorted.bam.sorted -fq "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".f2.fq -fq2 "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".f1.fq

# Add check to see whether both ref and personal alignments should be done, or just one (i.e. individuals not in UK10K)
../software/trim_galore_zip/trim_galore -stringency 5 -q 1 -o "$SAMPLE_ID" --phred33 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --paired "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".f1.fq "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".f2.fq
# Add trimming and adapter removal step.
perl ../software/prinseq-lite-0.20.4/prinseq-lite.pl -fastq "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".f1_val_1.fq -fastq2 "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".f2_val_2.fq -out_good "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID" -trim_tail_left 5 -trim_tail_right 5 -min_len 20

rm "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".f2_val_2.fq
rm "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".f1_val_1.fq
rm "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".f2.fq_trimming_report.txt
rm "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".f1.fq_trimming_report.txt
rm "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID"_1_singletons.fastq
rm "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID"_2_singletons.fastq
rm "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID"_prinseq_bad*
rm "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".f1.fq
rm "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".f2.fq
rm "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID"_sorted.bam.sorted # Ensure this is the name after replacing BAM_IDs with TwinUK IDs.
# Maternal genome generation (suffix arrays etc)
../software/STAR/bin/Linux_x86_64/STAR --runThreadN $THREAD_NO --runMode genomeGenerate --genomeDir "$TISSUE"/"$SAMPLE_ID"/maternal --genomeFastaFiles "$TISSUE"/"$SAMPLE_ID"/maternal/chr*_"$SAMPLE_ID"_maternal.fa

#Paternal genome generation (suffix arrays etc)
../software/STAR/bin/Linux_x86_64/STAR --runThreadN "$TISSUE"/$THREAD_NO --runMode genomeGenerate --genomeDir "$TISSUE"/"$SAMPLE_ID"/paternal --genomeFastaFiles "$TISSUE"/"$SAMPLE_ID"/paternal/chr*_"$SAMPLE_ID"_paternal.fa

#Align paternal
../software/STAR/bin/Linux_x86_64/STAR --runThreadN $THREAD_NO --runMode alignReads --readFilesIn "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID"_1.fastq "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID"_2.fastq --genomeDir "$TISSUE"/"$SAMPLE_ID"/paternal --outSAMstrandField intronMotif --outFilterMultimapNmax 30 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 15 --outMultimapperOrder Random --outSAMunmapped Within --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --sjdbOverhang 48 --outFilterMismatchNmax 6 --outSAMattributes NH nM NM MD HI --outSAMattrRGline  ID:"$SAMPLE_ID"_paternal PU:Illumina PL:Illumina LB:"$SAMPLE_ID"_paternal SM:"$SAMPLE_ID"_paternal CN:Seq_centre --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$SAMPLE_ID"_pat.
mv "$TISSUE"/"$SAMPLE_ID"_pat.Aligned.sortedByCoord.out.bam "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID"_pat.Aligned.sortedByCoord.out.bam
mv "$TISSUE"/"$SAMPLE_ID"_pat.Chimeric.out.junction "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID"_pat.Chimeric.out.junction
mv "$TISSUE"/"$SAMPLE_ID"_pat.Chimeric.out.sam "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID"_pat.Chimeric.out.sam
mv "$TISSUE"/"$SAMPLE_ID"_pat.Log.final.out "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID"_pat.Log.final.out
mv "$TISSUE"/"$SAMPLE_ID"_pat.SJ.out.tab "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID"_pat.SJ.out.tab
rm "$TISSUE"/"$SAMPLE_ID"_pat.Log.out
rm "$TISSUE"/"$SAMPLE_ID"_pat.Log.progress.out
#Align maternal
../software/STAR/bin/Linux_x86_64/STAR --runThreadN $THREAD_NO --runMode alignReads --readFilesIn "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID"_1.fastq "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID"_2.fastq --genomeDir "$TISSUE"/"$SAMPLE_ID"/maternal --outSAMstrandField intronMotif --outFilterMultimapNmax 30 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 15 --outMultimapperOrder Random --outSAMunmapped Within --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --sjdbOverhang 48 --outFilterMismatchNmax 6 --outSAMattributes NH nM NM MD HI --outSAMattrRGline  ID:"$SAMPLE_ID"_maternal PU:Illumina PL:Illumina LB:"$SAMPLE_ID"_maternal SM:"$SAMPLE_ID"_maternal CN:Seq_centre --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$SAMPLE_ID"_mat.
mv "$TISSUE"/"$SAMPLE_ID"_mat.Aligned.sortedByCoord.out.bam "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID"_mat.Aligned.sortedByCoord.out.bam
mv "$TISSUE"/"$SAMPLE_ID"_mat.Chimeric.out.junction "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID"_mat.Chimeric.out.junction
mv "$TISSUE"/"$SAMPLE_ID"_mat.Chimeric.out.sam "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID"_mat.Chimeric.out.sam
mv "$TISSUE"/"$SAMPLE_ID"_mat.Log.final.out "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID"_mat.Log.final.out
mv "$TISSUE"/"$SAMPLE_ID"_mat.SJ.out.tab "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID"_mat.SJ.out.tab
rm "$TISSUE"/"$SAMPLE_ID"_mat.Log.out
rm "$TISSUE"/"$SAMPLE_ID"_mat.Log.progress.out

# Add standard reference alignment too.

../software/STAR/bin/Linux_x86_64/STAR --runThreadN $THREAD_NO --runMode alignReads --readFilesIn "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID"_1.fastq "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID"_2.fastq --genomeDir hg19 --outSAMstrandField intronMotif --outFilterMultimapNmax 30 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 15 --outMultimapperOrder Random --outSAMunmapped Within --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --sjdbOverhang 48 --outFilterMismatchNmax 6 --outSAMattributes NH nM NM MD HI --outSAMattrRGline  ID:"$SAMPLE_ID"_maternal PU:Illumina PL:Illumina LB:"$SAMPLE_ID"_maternal SM:"$SAMPLE_ID"_maternal CN:Seq_centre --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$SAMPLE_ID"_ref.

mv "$TISSUE"/"$SAMPLE_ID"_ref.Aligned.sortedByCoord.out.bam "$TISSUE"/"$SAMPLE_ID"/reference/"$SAMPLE_ID"_ref.Aligned.sortedByCoord.out.bam
mv "$TISSUE"/"$SAMPLE_ID"_ref.Chimeric.out.junction "$TISSUE"/"$SAMPLE_ID"/reference/"$SAMPLE_ID"_ref.Chimeric.out.junction
mv "$TISSUE"/"$SAMPLE_ID"_ref.Chimeric.out.sam "$TISSUE"/"$SAMPLE_ID"/reference/"$SAMPLE_ID"_ref.Chimeric.out.sam
mv "$TISSUE"/"$SAMPLE_ID"_ref.Log.final.out "$TISSUE"/"$SAMPLE_ID"/reference/"$SAMPLE_ID"_ref.Log.final.out
mv "$TISSUE"/"$SAMPLE_ID"_ref.SJ.out.tab "$TISSUE"/"$SAMPLE_ID"/reference/"$SAMPLE_ID"_ref.SJ.out.tab
rm "$TISSUE"/"$SAMPLE_ID"_ref.Log.out
rm "$TISSUE"/"$SAMPLE_ID"_ref.Log.progress.out
########################################
#

../software/samtools-1.3.1/samtools view -b -F4 -q 30 "$TISSUE"/"$SAMPLE_ID"/reference/"$SAMPLE_ID"_ref.Aligned.sortedByCoord.out.bam -o "$TISSUE"/"$SAMPLE_ID"/reference/"$SAMPLE_ID".filtered.bam

## Produce consensus bams, in which the best read per haplotype is selected.
python ../software/PersonalGenomePipeline/pipeline.seperateBAMs.py "$SAMPLE_ID"
#
#for all variants without an rsid - assign them 'chr-pos'.
python ../software/PersonalGenomePipeline/nameVariants.py "$SAMPLE_ID" "$SAMPLE_ID".hets.phased.vcf.gz "$SAMPLE_ID".hets.phased.vcf.gz2

mv "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".hets.phased.vcf.gz2 "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".hets.phased.vcf.gz

# Sort the VCFs else picard/GATK throws a fit.
../software/vcftools/src/perl/vcf-sort -c "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".hets.phased.vcf.gz > "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.vcf

pigz --best -k "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.vcf

rm "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".hets.phased.vcf.gz
rm "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.vcf
#########################################
# run BASH script to concatenate reference genomes, and rename chain 1_maternal > 1

cp "$TISSUE"/"$SAMPLE_ID"/paternal.chain "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".paternal.edit.chain
cp "$TISSUE"/"$SAMPLE_ID"/maternal.chain "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".maternal.edit.chain
for i in {1..22}
do
	sed -i 's/'"$i"'_maternal/'"$i"'/g' "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".maternal.edit.chain
	sed -i 's/'"$i"'_paternal/'"$i"'/g' "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".paternal.edit.chain
	cat "$TISSUE"/"$SAMPLE_ID"/paternal/chr"$i"_"$SAMPLE_ID"_paternal.fa >> "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".all_paternal.fa
	cat "$TISSUE"/"$SAMPLE_ID"/maternal/chr"$i"_"$SAMPLE_ID"_maternal.fa >> "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".all_maternal.fa

done

sed -i 's/X_paternal/X/g' "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".paternal.edit.chain
#sed -i 's/Y_paternal/Y/g' "$SAMPLE_ID"/"$SAMPLE_ID".paternal.edit.chain
sed -i 's/M_paternal/MT/g' "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".paternal.edit.chain

sed -i 's/X_maternal/X/g' "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".maternal.edit.chain
#sed -i 's/Y_maternal/Y/g' "$SAMPLE_ID"/"$SAMPLE_ID".maternal.edit.chain
sed -i 's/M_maternal/MT/g' "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".maternal.edit.chain

cat "$TISSUE"/"$SAMPLE_ID"/paternal/chrX_"$SAMPLE_ID"_paternal.fa >> "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".all_paternal.fa
#cat "$SAMPLE_ID"/paternal/chrY_"$SAMPLE_ID"_paternal.fa >> "$SAMPLE_ID"/paternal/"$SAMPLE_ID".all_paternal.fa
cat "$TISSUE"/"$SAMPLE_ID"/paternal/chrM_"$SAMPLE_ID"_paternal.fa >> "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".all_paternal.fa

cat "$TISSUE"/"$SAMPLE_ID"/maternal/chrX_"$SAMPLE_ID"_maternal.fa >> "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".all_maternal.fa
#cat "$SAMPLE_ID"/maternal/chrY_"$SAMPLE_ID"_maternal.fa >> "$SAMPLE_ID"/maternal/"$SAMPLE_ID".all_maternal.fa
cat "$TISSUE"/"$SAMPLE_ID"/maternal/chrM_"$SAMPLE_ID"_maternal.fa >> "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".all_maternal.fa

#################################

# rename maternal and paternal fasta files with correct scaffold/chr name.
python ../software/PersonalGenomePipeline/renameFaChr.py "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".all_maternal.fa "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".maternal.renamed.fa "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".all_paternal.fa "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".paternal.renamed.fa

pigz -d "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.vcf.gz

# Modify vcf so all chromosomes are prefixed with 'chr'
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.vcf > "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.withChr.vcf

# At this point, the VCF, all.fasta and chain files are prefixed with 'chr' - so they are compatible.

mv "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.withChr.vcf "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.vcf

pigz --best -k "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.vcf

# Remap VCF to parental genome coordinates using modified CrossMap (fixedBugs)
python ../software/CrossMap-0.2.3/bin/CrossMap.py vcf "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".maternal.edit.chain "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.vcf "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".maternal.renamed.fa "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".maternal.vcf

python ../software/CrossMap-0.2.3/bin/CrossMap.py vcf "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".paternal.edit.chain "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.vcf "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".paternal.renamed.fa "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".paternal.vcf

# This script swaps the REF and ALT alleles according to whether it's the maternal or paternal haplotype.
python ../software/PersonalGenomePipeline/haplotypeVCFAlleles.py "$SAMPLE_ID" "$SAMPLE_ID".maternal.vcf "$SAMPLE_ID".paternal.vcf
#
mv "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".maternal.vcf.maternal2.vcf "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".maternal.vcf
mv "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".paternal.vcf.paternal2.vcf "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".paternal.vcf

java -jar /home/centos/scratch/software/picard-tools-2.4.1/picard.jar CreateSequenceDictionary R=""$TISSUE"/"$SAMPLE_ID"/paternal/$SAMPLE_ID".paternal.renamed.fa O="$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".paternal.renamed.dict
java -jar /home/centos/scratch/software/picard-tools-2.4.1/picard.jar CreateSequenceDictionary R=""$TISSUE"/"$SAMPLE_ID"/maternal/$SAMPLE_ID".maternal.renamed.fa O="$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".maternal.renamed.dict

java -jar /home/centos/scratch/software/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I=""$TISSUE"/$SAMPLE_ID"/maternal/"$SAMPLE_ID".consensus.mat.filtered.sorted.bam O=""$TISSUE"/$SAMPLE_ID"/maternal/"$SAMPLE_ID".consensus.mat.filtered.sorted.readGroup.bam  RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
java -jar /home/centos/scratch/software/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I=""$TISSUE"/$SAMPLE_ID"/paternal/"$SAMPLE_ID".consensus.pat.filtered.sorted.bam O=""$TISSUE"/$SAMPLE_ID"/paternal/"$SAMPLE_ID".consensus.pat.filtered.sorted.readGroup.bam  RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20

# both say maternal for chromosome because maternal was used as a template BAM when selecting best reads.
~/scratch/software/samtools-1.3.1/samtools view -h "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".consensus.mat.filtered.sorted.readGroup.bam  | sed -e 's/_maternal//g' >> "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".consensus.mat.filtered.sorted2.bam
mv "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".consensus.mat.filtered.sorted2.bam "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".consensus.mat.filtered.sorted.readGroup.bam

~/scratch/software/samtools-1.3.1/samtools view -h "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".consensus.pat.filtered.sorted.readGroup.bam  | sed -e 's/_maternal//g' >> "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".consensus.pat.filtered.sorted2.bam
mv "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".consensus.pat.filtered.sorted2.bam "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".consensus.pat.filtered.sorted.readGroup.bam

java -jar /home/centos/scratch/software/picard-tools-2.4.1/picard.jar ReorderSam I="$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".consensus.pat.filtered.sorted.readGroup.bam O="$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".consensus.pat.filtered.sorted.readGroup.bam2 R="$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".paternal.renamed.fa ALLOW_CONTIG_LENGTH_DISCORDANCE=true
java -jar /home/centos/scratch/software/picard-tools-2.4.1/picard.jar ReorderSam I="$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".consensus.mat.filtered.sorted.readGroup.bam O="$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".consensus.mat.filtered.sorted.readGroup.bam2 R="$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".maternal.renamed.fa ALLOW_CONTIG_LENGTH_DISCORDANCE=true

mv "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".consensus.pat.filtered.sorted.readGroup.bam2 "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".consensus.pat.filtered.sorted.readGroup.bam
mv "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".consensus.mat.filtered.sorted.readGroup.bam2 "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".consensus.mat.filtered.sorted.readGroup.bam

~/scratch/software/samtools-1.3.1/samtools index "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".consensus.mat.filtered.sorted.readGroup.bam
~/scratch/software/samtools-1.3.1/samtools index "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".consensus.pat.filtered.sorted.readGroup.bam

java -jar ~/scratch/software/GenomeAnalysisTK.jar -R "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".maternal.renamed.fa -T ASEReadCounter -o "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".ASE.mat.csv -I "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".consensus.mat.filtered.sorted.readGroup.bam -sites "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".maternal.vcf -dels -U ALLOW_N_CIGAR_READS -S SILENT
java -jar ~/scratch/software/GenomeAnalysisTK.jar -R "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".paternal.renamed.fa -T ASEReadCounter -o "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".ASE.pat.csv -I "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".consensus.pat.filtered.sorted.readGroup.bam -sites "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".paternal.vcf -dels -U ALLOW_N_CIGAR_READS -S SILENT

'''Currently, the ASE output is relative to the haplotype reference used (Maternal reference or paternal reference) 
   i.e. REF in ASE.mat.csv = maternal allele etc.To assess REF Bias (which we should not expect due to using Personal genomes), 
   and for easier interpretation, its better to have ASE in terms of universal REF / ALT. i.e. adding up all the REF/ALT counts 
   for each haplotype, consolidated into a single ASE file per individual.'''

python ../software/PersonalGenomePipeline/ASERefCord.py "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".ASE.mat.csv "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".ASE.pat.csv "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".hets.GATK.sorted.vcf "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".maternal.vcf "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".paternal.vcf "$TISSUE"/"$SAMPLE_ID".maternal.alleles.csv "$TISSUE"/"$SAMPLE_ID".paternal.alleles.csv "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".maternal.ref.csv "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".paternal.ref.csv
#################################

Rscript ../software/PersonalGenomePipeline/ASERefOut.R "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".maternal.ref.csv "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".paternal.ref.csv "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".ASE.csv

### END OF ASE PIPELINE ###

### START OF GENE-LEVEL COUNT GENERATION ###

../software/liftOver -gff gencode.v19.annotation.gtf "$TISSUE"/"$SAMPLE_ID"/maternal.chain "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".mat.gtf "$TISSUE"/"$SAMPLE_ID".not_lifted.txt -minMatch=0.00000001
../software/liftOver -gff gencode.v19.annotation.gtf "$TISSUE"/"$SAMPLE_ID"/paternal.chain "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".pat.gtf "$TISSUE"/"$SAMPLE_ID".not_lifted.txt -minMatch=0.00000001

# remove '_parental' string from the chromosome column in the liftedOver GTFs so it matches with the BAM files.
 sed -i 's/_maternal//g' "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".mat.gtf
 sed -i 's/_paternal//g' "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".pat.gtf

# Gene-level counts without multi-mapping - for each haplotype. 

# Haplotypic feature count (gene level counts)
../software/subread-1.5.0-p3-Linux-x86_64/bin/featureCounts -p -T 8 -a "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".mat.gtf -o "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".GeneCount_Mat.txt "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".consensus.mat.filtered.sorted.readGroup.bam
../software/subread-1.5.0-p3-Linux-x86_64/bin/featureCounts -p -T 8 -a "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".pat.gtf -o "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".GeneCount_Pat.txt "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".consensus.pat.filtered.sorted.readGroup.bam
../software/subread-1.5.0-p3-Linux-x86_64/bin/featureCounts -p -T 8 -a gencode.v19.annotation.gtf -o "$TISSUE"/"$SAMPLE_ID"/reference/"$SAMPLE_ID".GeneCount_Ref.txt "$TISSUE"/"$SAMPLE_ID"/reference/"$SAMPLE_ID".filtered.bam

Rscript ../software/PersonalGenomePipeline/AddHaploCounts.R "$SAMPLE_ID" "$TISSUE"/"$SAMPLE_ID"/maternal/"$SAMPLE_ID".GeneCount_Mat.txt "$TISSUE"/"$SAMPLE_ID"/paternal/"$SAMPLE_ID".GeneCount_Pat.txt "$TISSUE"/"$SAMPLE_ID"/"$SAMPLE_ID".GeneCount.Final.txt

