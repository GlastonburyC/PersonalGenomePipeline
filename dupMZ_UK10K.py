# Duplicate MZs and have a folder for each twin with their VCF.
import os,sys
# MZ_UK10K_Fat.toDuplicate.ids
with open(sys.argv[1],'r') as MZ:
	for line in (line.strip().split() for line in MZ):
		if line[0][-1:] == '2':
			twin = str(int(line[0])-1)
		else:
			twin = str(int(line[0])+1)
		os.system('cp -r '+twin+"/"+" "+line[0])
		os.system("mv "+line[0]+"/"+twin+".vcf.gz"+" "+line[0]+"/"+line[0]+".vcf.gz")
