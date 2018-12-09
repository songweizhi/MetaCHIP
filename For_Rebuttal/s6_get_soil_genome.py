import os

genome_meadata_file = '/srv/scratch/z5039045/MetaCHIP_rebuttal/nature10571_genome_metadata.txt'
genome_folder = '/srv/scratch/z5039045/MetaCHIP_rebuttal/downloaded_genomes_renamed'

Soil_all_genome_folder = '/srv/scratch/z5039045/MetaCHIP_rebuttal/Soil_all'
Soil_only_genome_folder = '/srv/scratch/z5039045/MetaCHIP_rebuttal/Soil_only'


# get genome list
Soil_all_taxon_list = set()
Soil_only_taxon_list = set()
for genome_metadata in open(genome_meadata_file):

    genome_metadata_split = genome_metadata.strip().split('\t')
    taxon_id = genome_metadata_split[1]
    niche = genome_metadata_split[2]

    if 'Soil' in niche:
        Soil_all_taxon_list.add(taxon_id)
        if niche == 'Non-Human, Soil':
            Soil_only_taxon_list.add(taxon_id)


print(len(Soil_all_taxon_list))
print(len(Soil_only_taxon_list))


# get Soil_all genomes
for soil_genome in Soil_all_taxon_list:
    pwd_soil_all_genome = '%s/taxon%s.fasta' % (genome_folder, soil_genome)
    if os.path.isfile(pwd_soil_all_genome) is True:
        cp_cmd = 'cp %s %s/' % (pwd_soil_all_genome, Soil_all_genome_folder)
        os.system(cp_cmd)


# get Soil_only genomes
for soil_only_genome in Soil_only_taxon_list:
    pwd_soil_only_genome = '%s/taxon%s.fasta' % (genome_folder, soil_only_genome)
    if os.path.isfile(pwd_soil_only_genome) is True:
        cp_cmd = 'cp %s %s/' % (pwd_soil_only_genome, Soil_only_genome_folder)
        os.system(cp_cmd)



