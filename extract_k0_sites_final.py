from Bio.Seq import Seq
import sys
from Bio import AlignIO
import os.path
from os.path import exists
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from itertools import combinations
import math
import pandas as pd
import re
import sys
from operator import itemgetter
from Bio.Seq import MutableSeq
from itertools import zip_longest


def list_generator(list):
	gapped_list = set()
	for element in list:
		gapped_list = gapped_list|gappify(element)
	return(gapped_list);


first_pos_k0 = set(['TAA'])
second_pos_k0 = set(['TTA','TTG','CTA','CGA','CTG','CGG','AGA','AGG'])
first_and_third_pos_k0 = set(['TGA'])
first_and_second_pos_k0 = set(['AAC','TTT','TCT','TAT','TGT','TTC','TCC','TAC','TGC','TCA','TCG','TAG','CTT','CCT','CAT','CGT','CTC','CCC','CAC','CGC','CCA','CAA','CCG','CAG','ATT','ACT','AAT','AGT','ATC','ACC','AGC','ATA','ACA','AAA','ACG','AAG','GTT','GCT','GAT','GGT','GTC','GCC','GAC','GGC','GTA','GCA','GAA','GGA','GTG','GCG','GAG','GGG'])
all_pos_k0 = set(['TGG','ATG'])


all_first= first_pos_k0 | all_pos_k0 | first_and_third_pos_k0 | first_and_second_pos_k0
all_second= second_pos_k0 | all_pos_k0 | first_and_second_pos_k0
all_first_third = first_and_third_pos_k0 | all_pos_k0

alignment = AlignIO.read(sys.argv[1],"fasta")
align = AlignIO.read(sys.argv[1],"fasta")
align_2 = AlignIO.read(sys.argv[1],"fasta")
align_3 = AlignIO.read(sys.argv[1],"fasta")
alignment_len = alignment.get_alignment_length()
#check if multipe of 3 or not
#append indices of codon positions that are k0
K0_id_matched_cols = set()
ids_to_ignore = set()
gaps = set()
#loop on the alignment ength with a step of three each time

for i in range(0, alignment_len, 3):
	new_aln = alignment[:,i:i + 3]
	new_aln_2 = alignment[:,i:i + 1]
	col = []
	nt = []
	for rec in new_aln:
		col.append(str(rec.seq))

	for j in range(0,len(col)):
		second_pos_only = all(elem.upper() in second_pos_k0 for elem in col)
		first_pos_only = all(elem.upper() in first_pos_k0 for elem in col)
		first_third_pos_only = all(elem.upper() in first_and_third_pos_k0 for elem in col)
		first_and_second_pos_only = all(elem.upper() in first_and_second_pos_k0 for elem in col)
		
		k0_all = all(elem.upper() in all_pos_k0 for elem in col)
		first_any = all(elem.upper() in all_first for elem in col)
		second_any = all(elem.upper() in all_second for elem in col)
		first_third_any = all(elem.upper() in all_first_third for elem in col)

		if k0_all:
			K0_id_matched_cols.add(i)
			K0_id_matched_cols.add(i+1)
			K0_id_matched_cols.add(i+2)	
		if first_third_any:
			K0_id_matched_cols.add(i)
			K0_id_matched_cols.add(i+2)
		if first_any: 
			if first_pos_only: 
				K0_id_matched_cols.add(i)
			else:
				K0_id_matched_cols.add(i)

		if second_any:

			if first_and_second_pos_only:
				K0_id_matched_cols.add(i)
				K0_id_matched_cols.add(i+1)
			else:
				K0_id_matched_cols.add(i+1)
gap_codons = sorted(gaps)

first_index_list = sorted(K0_id_matched_cols)

for record1, record2 in zip_longest(alignment,align_2) :
	record1.seq = MutableSeq(record1.seq)
	record2.seq = MutableSeq('')
	for item in first_index_list:
		record2.seq.append(record1.seq[item])
print(align_2.get_alignment_length())

AlignIO.write(align_2, "k0_no_within_codon_gap_"+ sys.argv[1], "fasta")
