# PersonalGenomePipeline

## Functionality

Given Two BAM files derived from mapping reads to Maternal and Paternal haplotypes, select the best read that maps to either, baed on mapping quality

# I have introduced a new SAM flag 'HT' which is prefixed with P/M (paternal/maternal) and followed by (i.e. _SRM)

# Rules
1. If reads map to the same position in either haplotype, select one at random (P_SRM or M_SRM)
2. If read maps to the same position, with a different MAPQ in one haplotype, select that highest MAPQ (P_SBM or M_SBM)
3. if read maps to a different position in either haplotype with the same MAPQ, choose one at random (P_DRM or M_DRM)
4. if read maps to a different position in either haplotypes with different MAPQ, choose the best MAPQ. (P_DBM or M_DBM)

*For example, P_SRM - paternally mapped read, read mapped to the same position in either haplotype, same MAPQ, chosen randomly
or M_DBM - maternally mapped read, read mapped to different positions in either haplotype, with different MAPQ, Maternal selected due to higher MAPQ score.
