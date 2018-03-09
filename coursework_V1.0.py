import Bio
import sys
import time
import os
from Bio import Entrez
from Bio import ExPASy
from Bio import SeqIO
from Bio import SwissProt
from Bio import Entrez
from Bio import SeqIO
from Bio import SeqRecord
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
from io import StringIO

class translate():
    def __init__(self,y):
        self.y=y
    def convert(self):
        return self.y.translate()

class overwrite():
    def __init__(self,filenames):
        self.filenames=filenames

    def merge(self):
        with open('final_seq.fasta', 'w') as outfile:
            for fname in self.filenames:
                with open(fname) as infile:
                    outfile.write(infile.read())

# main

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
protein_sequence=translate(y).convert()
rec1 = SeqRecord(id='Unkown', description='',seq=Seq(str(protein_sequence)))
myrecord = [rec1]
SeqIO.write(myrecord,'protein_sequence.fasta','fasta')

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


filenames = ['protein_sequence.fasta', 'MSA_als.fasta']

overwrite(filenames).merge()

cline = ClustalwCommandline("clustalw", infile="final_seq.fasta")
cline()

## build the tree
from Bio import Phylo
tree = Phylo.read("final_seq.dnd", "newick")
Phylo.draw_ascii(tree)

## delete useless files
os.remove('MSA_als.fasta')
os.remove('protein_sequence.fasta')
os.remove('final_seq.fasta')
os.remove('final_seq.dnd')
os.remove('final_seq.aln')
