import os

os.chdir('./blast_hits/')

con = open('blast_hits.txt')
use = con.read().split()

ids = list()
for i in range(2, len(use), 2):
	ids.append(use[i])

fin = list()
for i in range(0, len(ids)):
	fin.append(ids[i].replace('"', ''))

import Bio

from Bio import Entrez
from Bio import SeqIO

Entrez.email = 'cdotson@luc.edu'

hold1 = Entrez.efetch(db = 'nucleotide', id = fin, rettype = 'gb', retmode = 'text')
records = SeqIO.parse(hold1, "genbank") 
all_na = list()
for taxa in records:
	all_na.append(taxa.annotations['taxonomy'])

import csv

with open('ncbi_hits.csv', 'w') as f:
	write = csv.writer(f)
	write.writerows(all_na)
