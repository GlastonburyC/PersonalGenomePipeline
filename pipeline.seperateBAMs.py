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
import vcf

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
	for file in glob.glob(sys.argv[1]+'/'+"*.map"):
		chr=file.split('_')[0]
		Map_object[str(chr)]=MapParserInner(PARENT,chrom=file)
	#dictionary, key = chromosome, values are a dictionary of blocks.
	return Map_object


def translateMappedPosition(chr,cord,PARENT):
	# This function translates a reads mapping coordinate to that of the universal
	# reference genome. If the read maps within an insertion, it is assigned the 
	# last position unchanged in the reference (beginning of the insertion).
	ref_cord=[]
	if chr =='MT':
		ref_cord=cord
		return ref_cord
	if chr =='Y':
		ref_cord=cord
		return ref_cord
	if chr == '0':
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
def checkEditDistance(pedit,medit):
	if pedit == medit:
		return 0
	if pedit < medit:
		return 1
	else:
		return 2

def translateReadCord(mline,PARENT):
	# modify position
	if mline.pos == -1:
		chrom=0
	else:
		try:
			chrom = mline.reference_name
		except ValueError:
			chrom = 0
	mline.pos=translateMappedPosition(str(chrom),mline.pos,PARENT)
	# modify mate's position
	if mline.mpos == -1:
	    chrom=0
	else:
		try:
			chrom = mline.reference_name
		except ValueError:
			chrom = 0
	mline.mpos=translateMappedPosition(str(chrom),mline.mpos,PARENT)
	# calculate template length
	if mline.is_read1:
		mline.template_length=mline.mpos-mline.pos+49
	else:
		mline.template_length=mline.mpos-mline.pos-49
	return mline

def compareHapReads(mline,pline):
	flag=''
	if mline.qname == pline.qname:
		if checkEditDistance(pline.tags[2][1],mline.tags[2][1]) == 0:
			 y = random.random()
			 if y <0.5:
			 	mline.tags+=[('HT','random')]
			 	to_write=mline
			 	flag='M'
			 else:
			 	pline.tags+=[('HT','random')]
			 	to_write=pline
			 	flag='P'
		elif checkEditDistance(pline.tags[2][1],mline.tags[2][1]) == 1:
			pline.tags+=[('HT','best')]
			to_write=pline
			flag='P'
		else:
			mline.tags+=[('HT','best')]
			to_write=mline
			flag='M'
	return to_write,flag

def compareAndWrite(matr,patr,best_mat,best_pat):
	to_write,flag = compareHapReads(matr,patr)
	if flag=='M':
		best_mat.write(to_write)
	else:
		best_pat.write(to_write)


##### END OF FUNCTIONS #####

# Store maternal and paternal blocks relative to the reference.

maternal_map=MapParser(PARENT='M')
paternal_map=MapParser(PARENT='P')

os.system("/home/ubuntu/media/software/samtools-1.3.1/samtools sort -n "+sys.argv[1]+'/'+'maternal/'+sys.argv[1]+"_mat.Aligned.sortedByCoord.out.bam -o "+sys.argv[1]+'/'+'maternal/'+sys.argv[1]+"_mat.Aligned.sortedByCoord.out.sorted")
os.system("/home/ubuntu/media/software/samtools-1.3.1/samtools sort -n "+sys.argv[1]+'/'+'paternal/'+sys.argv[1]+"_pat.Aligned.sortedByCoord.out.bam -o "+sys.argv[1]+'/'+'paternal/'+sys.argv[1]+"_pat.Aligned.sortedByCoord.out.sorted")
os.system("mv "+sys.argv[1]+'/'+'paternal/'+sys.argv[1]+"_pat.Aligned.sortedByCoord.out.sorted "+sys.argv[1]+'/'+'paternal/'+sys.argv[1]+"_pat.Aligned.sortedByCoord.out.bam")
os.system("mv "+sys.argv[1]+'/'+'maternal/'+sys.argv[1]+"_mat.Aligned.sortedByCoord.out.sorted "+sys.argv[1]+'/'+'maternal/'+sys.argv[1]+"_mat.Aligned.sortedByCoord.out.bam")


# A new SAM flag is introduced HT = Haplotype. With this flag it's possible 
# to tell which haplotype a read mapped best to:
# P or M denote Paternal or Maternal Haplotype
# SBM = same mapping position, best MAPQ selected
# SRM = same mapping position, same MAPQ, randomly selected
# DBM = different mapping position, best MAPQ selected
# DRM = different mapping position, same MAPQ, randomly selected
# Output to Consensus.mat.bam or Consensus.pat.bam

mat = pysam.Samfile(sys.argv[1]+'/'+'maternal/'+sys.argv[1]+"_mat.Aligned.sortedByCoord.out.bam", 'rb')
mat_line_number = 0
for mline in mat.fetch(until_eof=True):
	mat_line_number += 1

pat = pysam.Samfile(sys.argv[1]+'/'+'paternal/'+sys.argv[1]+"_pat.Aligned.sortedByCoord.out.bam", 'rb')
pat_line_number = 0
for pline in pat.fetch(until_eof=True):
	pat_line_number += 1


#consensus = pysam.Samfile('consensus.bam','wb',template=mat)
mat = pysam.Samfile(sys.argv[1]+'/'+'maternal/'+sys.argv[1]+"_mat.Aligned.sortedByCoord.out.bam", 'rb')
pat = pysam.Samfile(sys.argv[1]+'/'+'paternal/'+sys.argv[1]+"_pat.Aligned.sortedByCoord.out.bam", 'rb')

best_mat = pysam.Samfile(sys.argv[1]+'/'+'maternal/'+sys.argv[1]+'.consensus.mat.bam','wb',template=mat)
best_pat = pysam.Samfile(sys.argv[1]+'/'+'paternal/'+sys.argv[1]+'.consensus.pat.bam','wb',template=mat)


matr = next(mat)
patr = next(pat)

i=1
j=1
# Iterate through both mat and pat simultaneously. Iterate through each read entry 
# (primary and secondary alignment) only extracting primary alignments
# Translate the pos and mpos of the read on the fly to reference coordinates.
# Write a consensus BAM file, consisting of the best read, given a haplotype alignment.
while (i < mat_line_number) and (j < pat_line_number):
	read_name = matr.qname
	patr_list = []
	matr_list = []
	# Get primary maternal reads
	while (matr.qname == read_name):
		if matr.is_secondary == False:
			matr_list.append(matr)
		if i < mat_line_number:
			matr = next(mat)
		else:
			break
		i+=1
	# Get primary paternal reads
	while (patr.qname == read_name):
		if patr.is_secondary == False:
			patr_list.append(patr)
		if j < pat_line_number:
			patr = next(pat)
		else:
			break
		j+=1
	if (len(patr_list) == 2) and (len(matr_list) == 2):
		matr1=matr_list[0]
		matr2=matr_list[1]
		patr1=patr_list[0]
		patr2=patr_list[1]
		compareAndWrite(matr1,patr1,best_mat,best_pat)
		compareAndWrite(matr2,patr2,best_mat,best_pat)
	elif (len(patr_list) == 1) and (len(matr_list) == 1):
		matr1=matr_list[0]
		patr1=patr_list[0]
		if (matr1.is_read1 and patr1.is_read1) or (matr1.is_read2 and patr1.is_read2):
			#matr1=translateReadCord(matr1,PARENT='M')
			#patr1=translateReadCord(patr1,PARENT='P')
			compareAndWrite(matr1,patr1,best_mat,best_pat)

best_pat.close()
best_mat.close()

ID=str(sys.argv[1])
vcf_reader = vcf.Reader(open(sys.argv[1]+'/'+sys.argv[1]+'.vcf.gz','r'))
vcf_writer = vcf.Writer(open(sys.argv[1]+'/'+sys.argv[1]+'.hets.phased.vcf','w'),vcf_reader)
for record in vcf_reader:
  if len(record.get_hets()) != 0:
  	if record.genotype(ID)['GT'][1]=='|':
  		if sorted(record.genotype(ID)['PL'])[1] >= 30:
  			vcf_writer.write_record(record)
vcf_writer.close()

# create het vcf index, else GATK throws a fit about Karyotypic ordering

# filter the BAM file

os.system('/media/shared_data/software/samtools-1.3.1/samtools view -b -F4 -q 30 '+sys.argv[1]+'/'+'maternal/'+sys.argv[1]+'.consensus.mat.bam -o '+sys.argv[1]+'/'+'maternal/'+sys.argv[1]+'.consensus.mat.filtered.bam')
os.system('/media/shared_data/software/samtools-1.3.1/samtools view -b -F4 -q 30 '+sys.argv[1]+'/'+'paternal/'+sys.argv[1]+'.consensus.pat.bam -o '+sys.argv[1]+'/'+'paternal/'+sys.argv[1]+'.consensus.pat.filtered.bam')


os.system('/media/shared_data/software/samtools-1.3.1/samtools sort '+sys.argv[1]+'/'+'paternal/'+sys.argv[1]+'.consensus.pat.filtered.bam -o '+sys.argv[1]+'/'+'paternal/'+sys.argv[1]+'.consensus.pat.filtered.sorted.bam')
os.system('/media/shared_data/software/samtools-1.3.1/samtools sort '+sys.argv[1]+'/'+'maternal/'+sys.argv[1]+'.consensus.mat.filtered.bam -o '+sys.argv[1]+'/'+'maternal/'+sys.argv[1]+'.consensus.mat.filtered.sorted.bam')

os.system('rm '+sys.argv[1]+'/'+'paternal/'+sys.argv[1]+'.consensus.pat.filtered.bam')
os.system('rm '+sys.argv[1]+'/'+'maternal/'+sys.argv[1]+'.consensus.mat.filtered.bam')

os.system('pigz --best -k '+sys.argv[1]+'/'+sys.argv[1]+'.hets.phased.vcf')
