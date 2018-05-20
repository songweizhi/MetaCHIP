
################################### run Get_clusters.py ##################################

python ~/PycharmProjects/MetaCHIP/Get_clusters.py -i Human_gut_bins -x fasta -p Human_gut
python ~/PycharmProjects/MetaCHIP/Get_clusters.py -i North_Sea_bins -x fasta -p North_Sea


####################################### run GTDB-Tk ######################################

gtdbtk classify_wf --genome_dir Human_gut_bins --out_dir Human_gut_bins_gtdbtk
gtdbtk classify_wf --genome_dir North_Sea_bins --out_dir North_Sea_bins_gtdbtk


################################### run get_homologues ###################################

get_homologues.pl -f 70 -t 3 -S 70 -E 1e-05 -C 70 -G -d Human_gut_input_genomes_gbk
get_homologues.pl -f 70 -t 3 -S 70 -E 1e-05 -C 70 -G -d North_Sea_input_genomes_gbk


################################### Best-match approach ##################################

python ~/PycharmProjects/MetaCHIP/Best-match.py -p Human_gut
python ~/PycharmProjects/MetaCHIP/Best-match.py -p North_Sea


################################## Phylogenetic approach #################################

python ~/PycharmProjects/MetaCHIP/Phylogenetic.py -p Human_gut
python ~/PycharmProjects/MetaCHIP/Phylogenetic.py -p North_Sea


##################################### COG annotation #####################################

# Human_gut bins, cd into Human_gut_faa_files folder produce by Get_clusters.py, then run
for each in *.faa; do python COG_wrapper.py -in $each -t P; done


# North_Sea bins, cd into North_Sea_faa_files folder produce by Get_clusters.py, then run
for each in *.faa; do python COG_wrapper.py -in $each -t P; done

# predicted HGTs from Human gut bins
python COG_wrapper.py -in Human_gut_HGT_candidates_Best-match_aa.fasta -t P
python COG_wrapper.py -in Human_gut_HGT_candidates_Phylogenetic_aa.fasta -t P

# predicted HGTs from North Sea bins
python COG_wrapper.py -in North_Sea_HGT_candidates_Best-match_aa.fasta -t P
python COG_wrapper.py -in North_Sea_HGT_candidates_Phylogenetic_aa.fasta -t P


############################ get the number of AR related COGs ###########################

# cd into COG annotation output folder, then run 
python get_AR_COGs.py
# path to type2cog.tab (line 13) and whog (line 14) files may need to be changed in the script. 
