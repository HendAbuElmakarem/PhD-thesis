#looping on different files to get the cds alignment
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from itertools import zip_longest
import itertools
import glob
import os
import shutil

dir_name= "/ceph/users/habu-elmakarem/hepatocystis_vinckeia/isolates_plus_pauls_data/Single_Copy_Orthologue_Sequences/mafft/cds_align/"
final_output = "/ceph/users/habu-elmakarem/hepatocystis_vinckeia/isolates_plus_pauls_data/Single_Copy_Orthologue_Sequences/mafft/cds_align/output/"
os.chdir(dir_name)
cds_files = sorted(glob.glob( 'cds*'))

mafft_files = sorted(glob.glob( 'mafft*'))

for (cds_file,mafft_file) in zip_longest(cds_files,mafft_files):
    #check that the sufix is the same in both files, i.e. we are dealing with the same orthogroup
    if cds_file.split('_')[0] == mafft_file.split("_")[1]:
        alignment = AlignIO.read(mafft_file,"fasta")
        alignment_len = alignment.get_alignment_length()
        pgon_cds_aln = SeqIO.parse(cds_file,"fasta")
        pgon_gene = SeqIO.parse(cds_file,"fasta")

        for cds in pgon_cds_aln:
            for record in alignment:
            #check that the species exists in both files
                if cds.id.split('_')[0] == record.id.split("_")[0]:
                    protein = []
                    gene  = []
                   if cds.seq.endswith('TGA') or cds.seq.endswith('TAG') or cds.seq.endswith('TAA'):
                       alignment_len = len(cds.seq)-3
                   else:
                       alignment_len = alignment.get_alignment_length()

                    for i in range(0,alignment_len,1):
                        protein.append(record[i])

                    for x in range(0,len(cds.seq),3):
                        gene.append(cds.seq[x:x+3])

                    for z, y in enumerate(protein):
                        if protein[z] == "-":
                            gene.insert(z,Seq("---"))
                    new_path = "/ceph/users/habu-elmakarem/hepatocystis_vinckeia/isolates_plus_pauls_data/Single_Copy_Orthologue_Sequences/mafft/cds_align/output/"
                    if os.path.isfile(os.path.join(new_path,cds_file.split('_')[1].split('.')[0]+"_cds_aln.fa")) == False:
                    #	print(os.path.join(path,cds_file.split('_')[1].split('.')[0]+"_cds_aln.fa"))
                        with open(cds_file.split('_')[0].split('.')[0]+"_cds_aln.fa","a") as out:
                                new_seq = Seq('')
                                for n in gene:
                                        new_seq +=n
                                out.write(">"+str(record.id) + "\n" + str(new_seq)+ "\n")
for file in glob.glob("*cds_aln*"):
    shutil.move(file,new_path)

