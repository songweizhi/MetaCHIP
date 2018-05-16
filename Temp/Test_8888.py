from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import SeqFeature as SF
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature


################ A: Make a SeqRecord ################

my_sequence = Seq("GATCGATCGATCGATCGATCGATCGATCGATC")
my_sequence_record = SeqRecord(my_sequence)
my_sequence_record.seq.alphabet = generic_dna

################ B: Make a SeqFeature ################

# define a FeatureLocation
my_start_pos = SF.ExactPosition(2)
my_end_pos = SF.ExactPosition(6)
my_feature_location = FeatureLocation(my_start_pos,my_end_pos, strand=-1)

# Define a feature type
my_feature_type = "CDS"
notes={'locus_tag':'bin014_00001', "translation":"MSFSSHCYMRPPSMLSV"}
# Create a SeqFeature
my_feature = SeqFeature(my_feature_location, type=my_feature_type, qualifiers=notes)

# Append Feature to SeqRecord
my_sequence_record.features.append(my_feature)

# print
print(my_sequence_record.format("gb"))

