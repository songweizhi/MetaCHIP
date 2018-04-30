import os
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


os.chdir('/Users/songweizhi/Desktop/refseq_test_20')

files = '*.fna'
file_list = [os.path.basename(file_name) for file_name in glob.glob(files)]

for each_genome in file_list:
    for each_seq in SeqIO.parse(each_genome, 'fasta'):
        if not ('plasmid' in each_seq.description):
            #print(each_seq.description[12:])

            strain_name = ''
            if each_seq.description[12:].endswith(', complete genome'):
                strain_name = each_seq.description[12:][:-17]
            elif each_seq.description[12:].endswith(' complete genome'):
                strain_name = each_seq.description[12:][:-16]
            elif each_seq.description[12:].endswith(', complete sequence'):
                strain_name = each_seq.description[12:][:-19]
            elif each_seq.description[12:].endswith(' complete sequence'):
                strain_name = each_seq.description[12:][:-18]
            else:
                strain_name = each_seq.description[12:]


            print(strain_name)


# remove chromosome from the end




        #print(each_seq.description)
        #print('%s\t%s' % (each_seq.id, each_seq.description))
        description_split = each_seq.description.split()
        sequence_status = ''
        if (each_seq.description[-15:] == 'complete genome') or (each_seq.description[-17:] == 'complete sequence'):
            sequence_status = 'completed'
        else:
            sequence_status = 'uncompleted'

        #print(sequence_status)




