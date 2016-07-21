# PersonalGenomePipeline

Executing the following will create SLURM sbatch submission scripts for each individual. It will create different scripts depending on whether the sample has both UK10K sequence data and RNA-seq, or just RNA-seq

``` ./workflow.sh all_samples UK10K_samples no_threads ```

## Functionality

Given Two BAM files derived from mapping reads to Maternal and Paternal haplotypes, select the best read that maps to either, based on edit distance.

For each read, simply choose the read which has the least number of mismatches, aligned to either haplotype. As both reads are the same, the alignment with the least mismatches should be considered the optimal choice.

0. Align to parental haplotypes using phased hets with GQ > 30. Also align to the reference genome for comparison.
1. Create paternal and maternal vcfs
2. run ASEReadCounter on both haplotypes (using each parental reference)- neccessary due to INDELs.
3. Correct bug in ASEReadCounter that outputs INDEL alleles incorrectly (i.e. Should be AAG/A, But the output is A/A)
4. Consolidate alleleCounts into REF/ALT according to universal reference vcf. 
5. Generate gene-level count files for maternal and paternal + reference genome. Merge the mat and pat haplotypes into one consensus gene    level count file.

##Output

1. Personalised genome Allele specific expression file.
2. Universal reference Allele specific expression file for comparison
2. Gene level counts using the personalised genomes 
3. Gene level counts using the universal reference genome

Each personalised sample weighs in at approximately 10Gb.
