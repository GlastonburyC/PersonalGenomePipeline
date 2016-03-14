import os
import sys
import pysam
from collections import defaultdict


pat = pysam.Samfile('Paternal.Aligned.sortedByCoord.out.bam', 'rb' )
mat = pysam.Samfile('Maternal.Aligned.sortedByCoord.out.bam', 'rb')

read_ids=[]

for rname in pat:
	print rname
	pat.mate(rname).qname

# Take a map file as input - create blocks for either haplotype to map read to
def MapParserInner(PARENT,chrom):
	par_start=[]
	par_end=[]
	ref_start=[]
	if PARENT == "M":
		par=2
	else:
		par=1
	with open(chrom,'r') as mapping_file:
		lines = mapping_file.readlines()
		for i in range(1,len(lines)):
			line=lines[i].split()
			if line[par] !='0':
				par_start.append(int(line[par]))
				ref_start.append(int(line[0]))
	for idx, val in enumerate(par_start):
		if idx == len(par_start)-1:
			par_end.append(int(par_start[idx]*2))
		else:
			par_end.append(int(par_start[idx+1])-1)
	MapInner = {'ref.start':ref_start,'par.start':par_start,'par.end':par_end}
	return MapInner

def MapParser(PARENT):
	Map_object={}
	for file in glob.glob("*.map"):
		chr=file.split('_')[0]
		Map_object[chr]=MapParserInner(PARENT,chrom=file)
	return Map_object


def translateMappedPosition(chr,cord,PARENT):
	ref_cord=[]
	if PARENT=="M":
		pat_map = maternal_map
	else:
		pat_map = paternal_map
	for i, block_start in enumerate(pat_map[chr]['par.start']):
		block_end=pat_map[chr]['par.end'][i]
		if block_start <= cord <= block_end:
			ref_cord = cord-pat_map[chr]['par.start'][i]+pat_map[chr]['ref.start'][i]
	return ref_cord

maternal_map=MapParser(PARENT='M')
paternal_map=MapParser(PARENT='P')


