import Bio
import sys
import time
import os
from Bio import Entrez
from Bio import ExPASy
from Bio import SeqIO
from Bio import SwissProt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio import SeqIO




## find the longest ORF
record = next(SeqIO.parse('/home/tianyu/CSC 8311/example.fasta', "fasta"))
x=record.seq
def ORF_finder(x):
    for i in range(0,len(x),1):
        startcodon = x[i:i+3]
        if startcodon == 'ATG':
            for j in range(i,len(x),3):
                stopcodons = x[j:j+3]
                sc=('TAA','TAG','TGA')
                if stopcodons in sc:
                    ORF=x[i:j+3]
                    return(ORF)
                    break

## translate
y=ORF_finder(x)
protein_sequence=y.translate()
f=open('protein_sequence.fasta','w+')
f.write('>UK_000000 Unknown protein [Unknown origanism]\n')
f.write(str(protein_sequence))
f.write('\n')
f.close

## blastp
result_handle = NCBIWWW.qblast("blastp", "nr", protein_sequence)
blast_record = NCBIXML.read(result_handle)
l=[]
for k in blast_record.alignments[:10]:
    blast_result=k.title
    l.append(k.title)


## find gi number
gi_number=[ele.split('|')[1] for ele in l]
print(gi_number)
## get fasta file from internet
Entrez.email = "t.cui3@newcastle.ac.uk"
create_file=open('MSA_als.fasta','w+')
create_file.close
for b in range(len(gi_number)):
    with Entrez.efetch(db="protein", rettype="fasta", retmode="text", id=gi_number[b]) as handle:
        search_results = handle.read()
        MSA_als=open('MSA_als.fasta','a+')
        MSA_als.write(search_results)
        MSA_als.write('\n')
        MSA_als.close

## merge two fasta files
files = ['protein_sequence.fasta', 'MSA_als.fasta']
with open('MSA_sequences.fasta', 'w') as outfile:
    for name in files:
        with open(name) as infile:
            outfile.write(infile.read())
            outfile.write('\n')

## MSA
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
from io import StringIO
from Bio import AlignIO
cline = ClustalwCommandline("clustalw", infile="MSA_sequences.fasta")
cline()

## build the tree
from Bio import Phylo
tree = Phylo.read("MSA_sequences.dnd", "newick")
Phylo.draw_ascii(tree)
