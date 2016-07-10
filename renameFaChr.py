'''Rename all chromosomes in fasta to reflect VCF
   using specified file, rename.txt'''
import sys

fasta= open(sys.argv[1])
newnames= open('/media/shared_data/software/PersonalGenomePipeline/rename.txt','r')
newfasta= open(sys.argv[2], 'w')

for line in fasta:
    if line.startswith('>'):
        newname= newnames.readline()
        newfasta.write(newname)
    else:
        newfasta.write(line)

fasta.close()
newnames.close()
newfasta.close()

fasta= open(sys.argv[3],'r')
newnames= open('/media/shared_data/software/PersonalGenomePipeline/rename.txt','r')
newfasta= open(sys.argv[4], 'w')

for line in fasta:
    if line.startswith('>'):
        newname= newnames.readline()
        newfasta.write(newname)
    else:
        newfasta.write(line)

fasta.close()
newnames.close()
newfasta.close()
