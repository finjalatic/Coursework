#Readme file
This is biologcal tool for building a phylogenetic tree of a unkonwn nucleotide sequence. It was written in python3.5.
To use this script, users need to download 'biopython' package. ((sudo) pip3 install biopython)
Some other tools need to be pre-downloaded as well.
  'clustalw'(sudo apt install clustalw) for MSA and phylogenetic tree
  'blast'(sudo apt install blast2) for blastp
Usrs can run this program in command line directly by typing 'python3 coursework.py yourfastafile.fasta number'. The first arguement should be your unknown sequence and the second one should be a number that how many branches showed in the tree.
For example, if you have a fasta file named 'example.fasta' and you want to see the tree with 10 branches, you can type in the terminal like this:
'python3 coursework.py example.fasta 10'
