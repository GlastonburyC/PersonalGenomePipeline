'''Rename all chromosomes in fasta to reflect VCF
   using specified file, rename.txt'''
import sys

fasta= open(sys.argv[0])
newnames= open('rename.txt')
newfasta= open(sys.argv[1], 'w')

for line in fasta:
    if line.startswith('>'):
        newname= newnames.readline()
        newfasta.write(newname)
    else:
        newfasta.write(line)

fasta.close()
newnames.close()
newfasta.close()

fasta= open(sys.argv[2])
newnames= open('rename.txt')
newfasta= open(sys.argv[3], 'w')

for line in fasta:
    if line.startswith('>'):
        newname= newnames.readline()
        newfasta.write(newname)
    else:
        newfasta.write(line)

fasta.close()
newnames.close()
newfasta.close()
