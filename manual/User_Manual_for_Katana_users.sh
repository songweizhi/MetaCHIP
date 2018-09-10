
# MetaCHIP manual (for Katana users)


############################## For installation (on Katana) ##############################

module load python/2.7.12
cd ~
mkdir mypythonenv
virtualenv --system-site-packages mypythonenv
. mypythonenv/bin/activate
pip install matplotlib
pip install ete3


##################################### For latter run #####################################

module load python/2.7.12
cd ~
. mypythonenv/bin/activate

module load prodigal/2.6.3
module load hmmer/3.1b2
module load mafft/7.310
module load fasttree/2.1.10
module load R/3.4.2


# Get_clusters.py
cd /srv/scratch/z5039045/MetaCHIP_demo
python /srv/scratch/z5039045/Softwares/MetaCHIP/Get_clusters.py -i input_genomes -x fna -p demo 


# Best-match.py
module load blast+/2.6.0
python /srv/scratch/z5039045/Softwares/MetaCHIP/Best-match.py -p demo -num_threads 9


# get_homologues.pl
cd demo_MetaCHIP_wd
/srv/scratch/z5039045/Softwares/gh/get_homologues-x86_64-20170918/get_homologues.pl -f 70 -t 3 -S 70 -E 1e-05 -C 70 -G -n 16 -d demo_gbk_files


# Phylogenetic.py
python /srv/scratch/z5039045/Softwares/MetaCHIP/Phylogenetic.py -p demo




