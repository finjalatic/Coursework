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

# get the amino acid sequence

## find the longest ORF
record = next(SeqIO.parse(sys.argv[1], "fasta"))
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
blast_number=sys.argv[2]
for k in blast_record.alignments[:int(blast_number)]:
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
class overwrite():
    def __init__(self,pro_seq):
        self.pro_seq=pro_seq

    def merge(self):
        f2=(open('MSA_als.fasta','a+'))
        f2.write('\n')
        f2.write(self.pro_seq)
        f2.close

f1=open('protein_sequence.fasta','r+')
pro_seq=str(f1.read())
f1.close

overwrite(pro_seq).merge()

## MSA
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
from io import StringIO
from Bio import AlignIO
cline = ClustalwCommandline("clustalw", infile="MSA_als.fasta")
cline()

## build the tree
from Bio import Phylo
tree = Phylo.read("MSA_als.dnd", "newick")
Phylo.draw_ascii(tree)

## delete useless files
os.remove('MSA_als.dnd',)
os.remove('MSA_als.aln')
os.remove('MSA_als.fasta')
os.remove('protein_sequence.fasta')
