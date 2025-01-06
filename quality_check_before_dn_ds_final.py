from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from itertools import zip_longest
import itertools
import glob
import os
from itertools import groupby

dir_name = os.getcwd()
os.chdir(dir_name)
cds_files = sorted(glob.glob('*cds_aln.fa'))
new_cds_files = sorted(glob.glob('*_new.fa'))
gaps_2 = []
cds_ids = []
discard_pile = set()
for cds_file in cds_files:
    pgon_cds_aln = SeqIO.parse(cds_file,"fasta")
    for cds in pgon_cds_aln:
        new_cds = []
        if cds.seq.endswith('TGA') or cds.seq.endswith('TAG') or cds.seq.endswith('TAA'):
           new_len = len(cds.seq)-3
        else:
           new_len = len(cds.seq)

        for n in range(0,new_len,3):
            new_cds.append(cds.seq[n:n+3])
        with open(cds_file.split('_')[0]+"_new.fa","a+") as c:
             with open(cds_file.split('_')[0]+"_new.fa") as out:
                  r = out.read()
                  if cds.id not in r:
                     with open(cds_file.split('_')[0]+"_new.fa","a") as out:
                          new_seq = Seq('')
                          for g in new_cds:
                              new_seq +=g
                          out.write(">"+str(cds.id) + "\n" + str(new_seq)+ "\n")

def all_equal(iterable):
    g = groupby(iterable)
    return next(g, True) and not next(g, False)

for cds_new in new_cds_files:
    new_aln = SeqIO.parse(cds_new,"fasta")
    lengths = []
    seq_len = {}
    seq_id = []
    gap_length = []
    gap_dict = {}
    gene_len = []
    for cd in new_aln:
        lengths.append(len(cd))
        seq_id.append(cd.id)
        for se in range(0,len(cd.seq),1):
            seq = []
            seq.append(cd.seq[se])
            for seq_list in seq:
                if seq_list == "-":
                   gap_length.append(seq_list)
                else:
                   gene_len.append(seq_list)

    for k ,v in zip(seq_id,lengths):
        seq_len[k] = v
    print(lengths)
    if all_equal(lengths) == False:
       print("sequences are unequal in length")
       {k:v for k, v in sorted(seq_len.items(), key=lambda kv: kv[1], reverse=True)}
       print("This sequence is longer than the others",list(seq_len.keys())[0])
       discard_pile.add(cds_new.split("_")[0])
    for key , val in seq_len.items():
        if val%3 !=0:
           print("sequence",key,"is not a multiple of three")
           discard_pile.add(cds_new.split("_")[0])

with open("discard_OGs.txt","a+") as f1:
    with open("discard_OGs.txt") as f2:
         ff = f2.read()
         for d in discard_pile:
             if d not in ff:
                with open("discard_OGs.txt","a") as f2:
                     f1.write(d+"\n")




