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
				par_start.append(int(line[par]))
				ref_start.append(int(line[0]))
	for idx, val in enumerate(par_start):
		# if last entry, set block to large number (should be end of chromosome)
		if idx == len(par_start)-1:
			par_end.append(int(par_start[idx]*2))
		else:
		# parental block end 	
			par_end.append(int(par_start[idx+1])-1)
	ref_end = np.array(ref_start) + (np.array(par_end)- np.array(par_start))
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
		ref_cord=cord =1
		return ref_cord
	ref_cord=[]
	if PARENT=="M":
		pat_map = maternal_map
	else:
		pat_map = paternal_map
	match=bisect.bisect(pat_map[chr]['par.blocks'], cord)
	if pat_map[chr]['ref.blocks'][match-1] == 0:
		ref_cord = pat_map[chr]['ref.blocks'][match-2]
	else:
		ref_cord = cord-pat_map[chr]['par.blocks'][match]+pat_map[chr]['ref.blocks'][match]
	return ref_cord

maternal_map=MapParser(PARENT='M')
paternal_map=MapParser(PARENT='P')


pat = pysam.Samfile('Paternal.Aligned.sortedByCoord.out.bam.sorted.bam', 'rb')
mat = pysam.Samfile('Maternal.Aligned.sortedByCoord.out.bam.sorted.bam', 'rb')

def TranslateAlignmentPos(PARENT):
	# for maternal and paternal alignments, store two dictionaries (key = mate id, values pos, mate.pos, qual, isize)
	# All positions are stored in reference coord space thanks to translateMappedPosition() to differentiate between reads
	# I have added whether the read is mate_1 or mate_2.
	read_name=[]
	pos=[]
	mpos=[]
	mapq=[]
	chromosome=[]
	isize=[]
	chrom=''
	if PARENT == "M":
		parent = mat
	else:
		parent = pat
	for read in parent.fetch(until_eof=True):
		if read.is_read1:
			qname=read.qname
			read_name.append(qname)
			if read.pos == -1:
				chromosome.append(0)
				pos.append(0)
				mpos.append(0)
				mapq.append(read.mapping_quality)
				isize.append(read.template_length)
			else:
				chrom=read.reference_name.split('_')[0]
				if chrom =="chrM":
					chromosome.append("chrM")
					pos.append(read.pos+1)
					mpos.append(read.mpos+1)
					mapq.append(read.mapping_quality)
					isize.append(read.mpos-read.pos+49)
				else:
					pos.append(translateMappedPosition(chrom,read.pos+1,PARENT=PARENT))
					mpos.append(translateMappedPosition(chrom,read.mpos+1,PARENT=PARENT))
					isize.append(read.mpos-read.pos+49)
					mapq.append(read.mapping_quality)
					chromosome.append(chrom)
		else:
			qname=read.qname
			read_name.append(qname)
			if read.pos == -1:
				chromosome.append(0)
				mapq.append(read.mapping_quality)
				isize.append(read.template_length)
				pos.append(read.pos)
				mpos.append(read.mpos)
			else:
				chrom=read.reference_name.split('_')[0]
				if chrom == "chrM":
					chromosome.append("chrM")
					pos.append(read.pos+1)
					mpos.append(read.mpos+1)
					mapq.append(read.mapping_quality)
					isize.append(read.mpos-read.pos-49)
				else:
					pos.append(translateMappedPosition(chrom,read.pos+1,PARENT=PARENT))
					mpos.append(translateMappedPosition(chrom,read.mpos+1,PARENT=PARENT))
					isize.append(read.mpos-read.pos-49)
					mapq.append(read.mapping_quality)
					chromosome.append(chrom)
	return read_name,chromosome,pos,mpos,mapq,isize


mat = pysam.Samfile('Maternal.Aligned.sortedByCoord.out.bam.sorted.bam', 'rb')
read_mat,chrom_mat,pos_mat,mpos_mat,mapq_mat,isize_mat = TranslateAlignmentPos(PARENT='M')

pat = pysam.Samfile('Paternal.Aligned.sortedByCoord.out.bam.sorted.bam', 'rb')
read_pat,chrom_pat,pos_pat,mpos_pat,mapq_pat,isize_pat = TranslateAlignmentPos(PARENT='P')


def samePosition(pqual,mqual):
	if pqual == mqual:
		return 0
	if pqual > mqual:
		return 1
	else:
		return 2


def diffPosition(pqual,mqual):
	if pqual == mqual:
		return 0
	if pqual > mqual:
		return 1
	else:
		return 2


def checkChromFormatting(chrom):
	try:
		chrom = mline.reference_name.split('_')[0]
	except ValueError:
		chrom = 0
	return chrom

consensus = pysam.Samfile('consensus.bam','wb',template=mat)
mat = pysam.Samfile('Maternal.Aligned.sortedByCoord.out.bam.sorted.bam', 'rb')
pat = pysam.Samfile('Paternal.Aligned.sortedByCoord.out.bam.sorted.bam', 'rb')
for mline in mat.fetch(until_eof=True):
	pline = next(pat)
	if mline.pos == -1:
		chrom=0
		mline.pos=translateMappedPosition(chrom,mline.pos+1,PARENT='M')
	else:
		try:
			chrom = mline.reference_name.split('_')[0]
		except ValueError:
			chrom = 0
		mline.pos=translateMappedPosition(chrom,mline.pos+1,PARENT='M')
	if mline.mpos == -1:
	    chrom=0
	    mline.mpos=translateMappedPosition(chrom,mline.mpos+1,PARENT='M')
	else:
		try:
			chrom = mline.reference_name.split('_')[0]
		except ValueError:
			chrom = 0
		mline.mpos=translateMappedPosition(chrom,mline.mpos+1,PARENT='M')
	if pline.pos == -1:
		chrom =0
		pline.pos=translateMappedPosition(chrom,pline.pos+1,PARENT='P')
	else: 
		try:
			chrom = pline.reference_name.split('_')[0]
		except ValueError:
			chrom = 0
		pline.pos=translateMappedPosition(chrom,pline.pos+1,PARENT='P')
	if pline.mpos == -1:
		chrom =0
		pline.mpos=translateMappedPosition(chrom,pline.mpos+1,PARENT='P')
	else:
		try:
			chrom = pline.reference_name.split('_')[0]
		except ValueError:
			chrom = 0
		pline.mpos=translateMappedPosition(chrom,pline.mpos+1,PARENT='P')
	if mread.is_read1:
		mread.template_length=mline.mpos-mline.pos+49
		pread.template_length=pline.mpos-pline.pos+49
	else:
		mread.template_length=mline.mpos-mline.pos-49
		pread.template_length=pline.mpos-pline.pos-49
	if mline.qname == pline.qname:
		if mline.pos == pline.pos:
			if samePosition(mline.mapping_quality,pline.mapping_quality) == 0:
				 y = random.random()
				 if y <0.5:
				 	consensus.write(mline)
				 else:
				 	consensus.write(pline)
			elif samePosition(mline.mapping_quality,pline.mapping_quality) == 1:
				consensus.write(pline)
			else:
				consensus.write(mline)
		if mline.pos != pline.pos:
			if diffPosition(pline.mapping_quality,mline.mapping_quality) == 0:
				y = random.random()
				if y <0.5:
					consensus.write(mline)
				else:
					consensus.write(pline)
			elif diffPosition(pline.mapping_quality,mline.mapping_quality) ==1:
				consensus.write(pline)
			else:
				consensus.write(mline)




mat = pysam.Samfile('Maternal.Aligned.sortedByCoord.out.bam.sorted.bam', 'rb')
pat = pysam.Samfile('Paternal.Aligned.sortedByCoord.out.bam.sorted.bam', 'rb')
for mline in mat.fetch(until_eof=True):
	pline = next(pat)
	if pline.qname == mline.qname:
		count = count+1
		print 'true '+str(count)
	else:
		print mline
		break
