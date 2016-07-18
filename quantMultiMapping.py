#### Add a column that counts the number of reads that map to more than 1 location

# quanyMultiMapping.py '$SAMPLE_ID'/maternal/'$SAMPLE_ID'.maternal.ref.csv '$SAMPLE_ID'/paternal/'$SAMPLE_ID'.paternal.ref.csv ASE.mat.ref.multi.txt ASE.pat.ref.multi.txt
import pysam,sys,os

consolidated_ase=open(sys.argv[1],'r')
header=next(consolidated_ase)
#1152_mat.Aligned.sortedByCoord.out.bam
mat = pysam.Samfile(sys.argv[5]+'/maternal/'+sys.argv[5]+'_mat.Aligned.sortedByCoord.out.bam', 'rb')
multi_out=open(sys.argv[3],'w')

multi_out.write(header.split('\n')[0]+"\t"+"MultiMaternal"+'\t'+"Maternal_totalReads"+'\n')
NH=0
total=0
for snp in consolidated_ase:
	line=snp.split('\t')
	chr=line[0]+"_maternal"
	for entry in mat.fetch(chr,int(line[1])-1,int(line[1])):
		if entry.tags[0][1] > 1:
			NH+=1
			total+=1
		else:
			total+=1
	multi_out.write("\t".join(line).split('\n')[0]+'\t'+str(NH)+'\t'+str(total)+"\n")
	NH=0
	total=0


multi_out.close()

consolidated_ase=open(sys.argv[2],'r')
multi_out=open(sys.argv[4],'w')
#1152_pat.Aligned.sortedByCoord.out.bam
pat = pysam.Samfile(sys.argv[5]+'/paternal/'+sys.argv[5]+'_pat.Aligned.sortedByCoord.out.bam', 'rb')

header=next(consolidated_ase)

multi_out.write(header.split('\n')[0]+"\t"+"MultiPaternal"+"\t"+"Paternal_totalReads"+'\n')
NH=0
total=0
for snp in consolidated_ase:
	line=snp.split('\t')
	chr=line[0]+"_paternal"
	for entry in pat.fetch(chr,int(line[1])-1,int(line[1])):
		if entry.tags[0][1] > 1:
			NH+=1
			total+=1
		else:
			total+=1
	multi_out.write("\t".join(line).split('\n')[0]+'\t'+str(NH)+'\t'+str(total)+"\n")
	NH=0
	total=0


multi_out.close()
