import os


# extract path to this config file
pwd_config_file = os.path.realpath(__file__)
config_file_path = '/'.join(pwd_config_file.split('/')[:-1])


# specify full path to corresponding executables at the right side of colon
config_dict = {'prodigal'        : 'prodigal',
               'hmmsearch'       : 'hmmsearch',
               'hmmfetch'        : 'hmmfetch',
               'hmmalign'        : 'hmmalign',
               'hmmstat'         : 'hmmstat',
               'mafft'           : 'mafft',
               'blastp'          : 'blastp',
               'blastn'          : 'blastn',
               'makeblastdb'     : 'makeblastdb',
               'fasttree'        : 'FastTree',
               'ranger_mac'      : '%s/Ranger-DTL-Dated.mac'   % config_file_path,  # do not edit this line
               'ranger_linux'    : '%s/Ranger-DTL-Dated.linux' % config_file_path,  # do not edit this line
               'path_to_hmm'     : '%s/MetaCHIP_phylo.hmm'     % config_file_path,  # do not edit this line
               'circos_HGT_R'    : '%s/MetaCHIP_circos_HGT.R'  % config_file_path   # do not edit this line
               }

