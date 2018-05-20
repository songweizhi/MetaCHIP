import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


'''
module load python/3.5.2
module load perl/5.20.1
module load blast+/2.6.0

python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in input.fasta -t N
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in input.fasta -t P

'''


def dna2aa(dna_file, aa_file):
    query_aa_handle = open(aa_file, 'w')
    for each in SeqIO.parse(dna_file, 'fasta'):
        each_aa = each.seq.translate()
        each_aa_record = SeqRecord(each_aa)
        each_aa_record.id = each.id
        each_aa_record.description = each.description
        SeqIO.write(each_aa_record, query_aa_handle, 'fasta')
    query_aa_handle.close()


parser = argparse.ArgumentParser()

parser.add_argument('-in',
                    required=True,
                    help='path to input sequences (in multi-fasta format)')

parser.add_argument('-t',
                    required=True,
                    help='The type of input sequences, "N" for "nucleotide", "P" for "protein"')

parser.add_argument('-rpsblast',
                    required=False,
                    default='rpsblast',
                    help='path to rpsblast executable')

parser.add_argument('-cdd2cog',
                    required=False,
                    default='/srv/scratch/z5039045/Scripts/cdd2cog.pl',
                    help='path to cdd2cog.pl file')

parser.add_argument('-db',
                    required=False,
                    default='/srv/scratch/z5039045/COG_annotation_db/Cog',
                    help='path to db file')

parser.add_argument('-fun',
                    required=False,
                    default='/srv/scratch/z5039045/COG_annotation_db/fun.txt',
                    help='path to fun.txt file')

parser.add_argument('-cddid',
                    required=False,
                    default='/srv/scratch/z5039045/COG_annotation_db/cddid.tbl',
                    help='path to cddid.tbl file')

parser.add_argument('-whog',
                    required=False,
                    default='/srv/scratch/z5039045/COG_annotation_db/whog',
                    help='path to whog file')

args = vars(parser.parse_args())
input_seq = args['in']
sequence_type = args['t']
pwd_rpsblast = args['rpsblast']
pwd_cdd2cog = args['cdd2cog']
pwd_db = args['db']
pwd_fun = args['fun']
pwd_cddid = args['cddid']
pwd_whog = args['whog']


# check whether file exist
unfound_inputs = []
for each_input in [pwd_cdd2cog, pwd_fun, pwd_cddid, pwd_whog]:
    if (not os.path.isfile(each_input)) and (not os.path.isdir(each_input)):
        unfound_inputs.append(each_input)
if len(unfound_inputs) > 0:
    for each_unfound in unfound_inputs:
        print('%s not found' % each_unfound)
    exit()


input_seq_no_ext, input_seq_ext = os.path.splitext(input_seq)
rpsblast_output = '%s_COG.tab' % input_seq_no_ext
output_folder = '%s_COG_results' % input_seq_no_ext


input_seq_aa = ''
if sequence_type == 'N':
    input_seq_aa = '%s_aa.fasta' % input_seq_no_ext
    dna2aa(input_seq, input_seq_aa)
if sequence_type == 'P':
    input_seq_aa = input_seq


print('Input file found: %s' % input_seq)
print('Running rpsblast ...')
os.system('%s -query %s -db %s -out %s -evalue 1e-2 -outfmt 6' % (pwd_rpsblast, input_seq_aa, pwd_db, rpsblast_output))
print('Running cdd2cog.pl ...')
os.system('perl %s -r %s -c %s -f %s -w %s' % (pwd_cdd2cog, rpsblast_output, pwd_cddid, pwd_fun, pwd_whog))
os.system('mv results %s' % (output_folder))

# remove temporary files
os.system('rm %s' % rpsblast_output)
if sequence_type == 'N':
    os.system('rm %s' % input_seq_aa)

print('All done, annotation results exported to %s' % output_folder)

