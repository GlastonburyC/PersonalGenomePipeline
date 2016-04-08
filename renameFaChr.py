'''Rename all chromosomes in fasta to reflect VCF
   using specified file, rename.txt'''

fasta= open('all_paternal.fa')
newnames= open('rename.txt')
newfasta= open('paternal.renamed.fa', 'w')

for line in fasta:
    if line.startswith('>'):
        newname= newnames.readline()
        newfasta.write(newname)
    else:
        newfasta.write(line)

fasta.close()
newnames.close()
newfasta.close()

fasta= open('all_maternal.fa')
newnames= open('rename.txt')
newfasta= open('maternal.renamed.fa', 'w')

for line in fasta:
    if line.startswith('>'):
        newname= newnames.readline()
        newfasta.write(newname)
    else:
        newfasta.write(line)

fasta.close()
newnames.close()
newfasta.close()
