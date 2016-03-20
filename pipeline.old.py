##### C. Glastonbury Haplotype coordinates to reference space BAM file creation ######
##### Requires latest version of PySam - Cluster currently has old version - this doesn't work! #####
import os
import sys
import pysam
from collections import defaultdict
import numpy as np
import glob
from itertools import chain
import bisect
import random

#### temporary function to test timing ####
if __name__=='__main__':
    from timeit import Timer
    t = Timer(lambda: translateMappedPosition('chr1',28592,PARENT='M'))
    print t.timeit(number=1)

# Take a map file as input - create blocks for either haplotype to map read to
def MapParserInner(PARENT,chrom):
	par_start=[]
	par_end=[]
	ref_start=[]
	ref_end=[]
	if PARENT == "M":
		par=2
	else:
		par=1
		# open chromosomal map file
	with open(chrom,'r') as mapping_file:
		lines = mapping_file.readlines()
		for i in range(1,len(lines)):
			line=lines[i].split()
			if line[par] !='0':
				# parental block start
				par_start.append(int(line[par])-1)
				ref_start.append(int(line[0])-1)
	for idx, val in enumerate(par_start):
		# if last entry, set block to large number (should be end of chromosome)
		if idx == len(par_start)-1:
			par_end.append(int(par_start[idx]*2))
		else:
		# parental block end 	
			par_end.append(int(par_start[idx+1])-1)
	for index,value in enumerate(ref_start):
		if ref_start[index]== -1:
			ref_end.append(-1)
		else:
			ref_end.append(ref_start[index] + (par_end[index] - par_start[index]))
	ref_blocks = zip(ref_start,ref_end)
	ref_blocks = list(chain(*ref_blocks))
	par_blocks = zip(par_start, par_end)
	par_blocks = list(chain(*par_blocks))
	# a dictionary of three lists reference start and end position, parental start and end position.
	MapInner = {'ref.blocks':ref_blocks,'par.blocks':par_blocks}
	return MapInner

def MapParser(PARENT):
	# iterate through all chromosomes 1-22, X, M
	Map_object={}
	for file in glob.glob("*.map"):
		chr=file.split('_')[0]
		Map_object[chr]=MapParserInner(PARENT,chrom=file)
	#dictionary, key = chromosome, values are a dictionary of blocks.
	return Map_object


def translateMappedPosition(chr,cord,PARENT):
	# This function translates a reads mapping coordinate to that of the universal
	# reference genome. If the read maps within an insertion, it is assigned the 
	# last position unchanged in the reference (beginning of the insertion).
	ref_cord=[]
	if chr =='chrM':
		ref_cord=cord
		return ref_cord
	if chr == 0:
		ref_cord = -1
		return ref_cord
	if cord == -1:
		ref_cord = -1
		return ref_cord
	if PARENT=="M":
		pat_map = maternal_map
	else:
		pat_map = paternal_map
	##
	match=bisect.bisect_left(pat_map[chr]['par.blocks'], cord)
	## even match denotes block start
	## odd match denotes block end
	## if cord is inside block, match will return end (odd number), so we get turn it into start
	if match % 2 != 0:
		match = match - 1
	ref_start = pat_map[chr]['ref.blocks'][match]
	par_start = pat_map[chr]['par.blocks'][match]
	if ref_start == -1:
		## map read to last observed reference position
		ref_cord = -1
		i=1
		while ref_cord == -1:
			ref_cord = pat_map[chr]['ref.blocks'][match-i]
			i+=2
	else:
		ref_cord = cord - par_start + ref_start
	return ref_cord

# reads in both haplotypes map to same position
def checkMAPQ(pqual,mqual):
	if pqual == mqual:
		return 0
	if pqual > mqual:
		return 1
	else:
		return 2


# keep only primary aligned reads
def PrimaryAlignedReads(bam_object):
	primary_reads=[]
	qname_out=[]
	for read in bam_object.fetch(until_eof=True):
		if read.is_secondary:
			pass
		else:
			if read.is_read1:
				qname=read.qname
				qname+='_1'
				read.qname=qname
				primary_reads.append(read)
				qname_out.append(qname)
			else:
				qname=read.qname
				qname+='_2'
				read.qname=qname
				primary_reads.append(read)
				qname_out.append(qname)
	return primary_reads,qname_out


def keepPrimaryAlignments(par_primary,primary):
	par_primary_out=[]
	for line in par_primary:
		if line.qname in primary:
			par_primary_out.append(line)
		else:
			pass
	return par_primary_out

##### END OF FUNCTIONS #####

# Store maternal and paternal blocks relative to the reference.
maternal_map=MapParser(PARENT='M')
paternal_map=MapParser(PARENT='P')

# Remove non-primary alignments
mat = pysam.Samfile('Maternal.Aligned.sortedByCoord.out.bam', 'rb')
mat_primary, mat_qname = PrimaryAlignedReads(mat)
pat = pysam.Samfile('Paternal.Aligned.sortedByCoord.out.bam', 'rb')
pat_primary, pat_qname = PrimaryAlignedReads(pat)

# It's possible one haplotype has more reads than the other after removing secondary alignments
# this is because split reads can have 1 or more primary alignments
# this difference is tiny, for this sample, the haplotype BAMs differed by 3 reads!
primary = set(mat_qname).intersection(pat_qname)

# Paternal primary alignments
pat_primary_out = keepPrimaryAlignments(pat_primary,primary)

# Maternal primary alignments
mat_primary_out = keepPrimaryAlignments(mat_primary,primary)

# write to BAM file
mat_primary_bam = pysam.Samfile('mat_primary.bam','wb',template=mat)
for line in mat_primary_out:
	mat_primary_bam.write(line)

# sort BAM file
mat_primary_bam.close()
os.system("samtools sort -n mat_primary.bam mat_primary_sorted")
os.system("rm mat_primary.bam")
mat_primary_out=None
# Write to BAM file
pat_primary_bam = pysam.Samfile('pat_primary.bam','wb',template=pat)
for line in pat_primary_out:
	pat_primary_bam.write(line)

pat_primary_bam.close()
pat_primary_out=None
# Sort BAM file

os.system("samtools sort -n pat_primary.bam pat_primary_sorted")
os.system("rm pat_primary.bam")


# This can be tidied up significantly, this is where the comparisons take place of reads mapping
# across two haplotypes. 
# A new SAM flag is introduced HT = Haplotype. With this flag it's possible 
# to tell which haplotype a read mapped best to:
# P or M denote Paternal or Maternal Haplotype
# SBM = same mapping position, best MAPQ selected
# SRM = same mapping position, same MAPQ, randomly selected
# DBM = different mapping position, best MAPQ selected
# DRM = different mapping position, same MAPQ, randomly selected
# Output to Consensus.bam
consensus = pysam.Samfile('consensus.bam','wb',template=mat)
mat = pysam.Samfile('mat_primary_sorted.bam', 'rb')
pat = pysam.Samfile('pat_primary_sorted.bam', 'rb')
for mline in mat.fetch(until_eof=True):
	pline = next(pat)
	if mline.pos == -1:
		chrom=0
		mline.pos=translateMappedPosition(chrom,mline.pos,PARENT='M')
	else:
		try:
			chrom = mline.reference_name.split('_')[0]
		except ValueError:
			chrom = 0
		mline.pos=translateMappedPosition(chrom,mline.pos,PARENT='M')
	if mline.mpos == -1:
	    chrom=0
	    mline.mpos=translateMappedPosition(chrom,mline.mpos,PARENT='M')
	else:
		try:
			chrom = mline.reference_name.split('_')[0]
		except ValueError:
			chrom = 0
		mline.mpos=translateMappedPosition(chrom,mline.mpos,PARENT='M')
	if pline.pos == -1:
		chrom =0
		pline.pos=translateMappedPosition(chrom,pline.pos,PARENT='P')
	else: 
		try:
			chrom = pline.reference_name.split('_')[0]
		except ValueError:
			chrom = 0
		pline.pos=translateMappedPosition(chrom,pline.pos,PARENT='P')
	if pline.mpos == -1:
		chrom =0
		pline.mpos=translateMappedPosition(chrom,pline.mpos,PARENT='P')
	else:
		try:
			chrom = pline.reference_name.split('_')[0]
		except ValueError:
			chrom = 0
		pline.mpos=translateMappedPosition(chrom,pline.mpos,PARENT='P')
	if mline.is_read1:
		mline.template_length=mline.mpos-mline.pos+49
		pline.template_length=pline.mpos-pline.pos+49
	else:
		mline.template_length=mline.mpos-mline.pos-49
		pline.template_length=pline.mpos-pline.pos-49
	if mline.qname == pline.qname:
		if mline.pos == pline.pos:
			if checkMAPQ(pline.mapping_quality,mline.mapping_quality) == 0:
				 y = random.random()
				 if y <0.5:
				 	mline.tags+=[('HT','M_SRM')]
				 	consensus.write(mline)
				 	print "Mat equal"
				 else:
				 	pline.tags+=[('HT','P_SRM')]
				 	print "Pat equal"
				 	consensus.write(pline)
			elif checkMAPQ(pline.mapping_quality,mline.mapping_quality) == 1:
				print "Pat equal greater"
				pline.tags+=[('HT','P_SBM')]
				consensus.write(pline)
			else:
				mline.tags+=[('HT','M_SBM')]
				print "Mat equal greater"
				consensus.write(mline)
		if mline.pos != pline.pos:
			if checkMAPQ(pline.mapping_quality,mline.mapping_quality) == 0:
				y = random.random()
				if y <0.5:
					print "Mat diff equal"
					mline.tags+=[('HT','M_DRM')]
					consensus.write(mline)
				else:
					print "Pat diff equal"
					pline.tags+=[('HT','P_DRM')]
					consensus.write(pline)
			elif checkMAPQ(pline.mapping_quality,mline.mapping_quality) ==1:
				print "Pat diff greater"
				pline.tags+=[('HT','P_DBM')]
				consensus.write(pline)
			else:
				mline.tags+=[('HT','M_DBM')]
				print "Mat diff greater"
				consensus.write(mline)
# File close to print E0F byte.
consensus.close()
os.system("samtools sort consensus.bam consensus.sorted")
os.system("rm consensus.bam")