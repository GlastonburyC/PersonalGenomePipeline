import os
import sys
import pysam
from collections import defaultdict
import numpy as np
import glob

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
	# a dictionary of three lists reference start and end position, parental start and end position.
	MapInner = {'ref.start':ref_start,'ref.end':ref_end,'par.start':par_start,'par.end':par_end}
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
	# last position unchanged in the reference (beginning of insertion).
	ref_cord=[]
	if PARENT=="M":
		pat_map = maternal_map
	else:
		pat_map = paternal_map
	for i, block_start in enumerate(pat_map[chr]['par.start']):
		block_end=pat_map[chr]['par.end'][i]
		if block_start <= cord <= block_end:
			# When considering indels - translate the position to the ref coordinate before indel
			if pat_map[chr]['ref.start'][i] == 0:
				ref_cord = pat_map[chr]['ref.end'][i-1]	
			else:
				ref_cord = cord-pat_map[chr]['par.start'][i]+pat_map[chr]['ref.start'][i]
	return ref_cord

maternal_map=MapParser(PARENT='M')
paternal_map=MapParser(PARENT='P')


pat = pysam.Samfile('paternal/Paternal.Aligned.sortedByCoord.out.bam', 'rb')
mat = pysam.Samfile('maternal/Maternal.Aligned.sortedByCoord.out.bam', 'rb')


def TranslateAlignmentPos(PARENT):
	parental2ref={}
	chr=''
	if PARENT == "M":
		parent = mat
	else:
		parent = pat
	for read in parent:
		if read.is_read1:
			qname=read.qname+'_1'
			if read.pos == -1:
				chr=0
				parental2ref[qname] = 0, read.mapping_quality	
			else:
				chr=read.reference_name.split('_')[0]
				parental2ref[qname] = chr, translateMappedPosition(chr,read.pos+1,PARENT=PARENT), read.mapping_quality, read.template_length
		else:
			qname=read.qname+'_2'
			if read.pos == -1:
				chr=0
				parental2ref[qname] = 0, read.mapping_quality
			else:
				chr=read.reference_name.split('_')[0]
				parental2ref[qname] = chr, translateMappedPosition(chr,read.pos+1,PARENT=PARENT), read.mapping_quality, read.template_length
	return parental2ref


mat2ref = TranslateAlignmentPos(PARENT='M')
pat2ref = TranslateAlignmentPos(PARENT='P')



