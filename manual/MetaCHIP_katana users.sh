
#################################### install MetaCHIP ####################################

module load python/2.7.15
cd ~
mkdir mypythonenv_metachip
virtualenv --system-site-packages mypythonenv_metachip
. mypythonenv_metachip/bin/activate
pip install MetaCHIP

# and make you have the following R packages installed
# optparse: command line option parser in R.
# ape: package for analyses of phylogenetics and evolution in R.
# circlize: package for circular visualization.

###################################### run MetaCHIP ######################################

#!/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l vmem=60gb
#PBS -l walltime=11:59:00
#PBS -j oe
#PBS -m ae

module unload python/3.5.2
module load python/2.7.15
cd ~
. mypythonenv_metachip/bin/activate
module load hmmer/3.1b2
module load mafft/7.310
module load fasttree/2.1.10
module load R/3.4.2
module load blast+/2.6.0
module load prodigal/2.6.3
export PATH=/srv/scratch/z5039045/Softwares/Ranger-DTL2/Linux/CorePrograms:$PATH
export PATH=/srv/scratch/z5039045/Softwares/Ranger-DTL2/Linux/SupplementaryPrograms:$PATH

cd /srv/scratch/z5039045/MetaCHIP/NorthSea_67GoodBins

# look at HGTs at Class level
MetaCHIP PI -i NorthSea_67GoodBins -x fasta -taxon NorthSea_67GoodBins.bac120.summary.tsv -p NorthSea_67GoodBins -r c -t 12 -qsub 

# run the following two steps after all submitted jobs were finished, 
# you can check the running status of your submitted jobs with qstat.
MetaCHIP BM -p NorthSea_67GoodBins -r c -t 12
MetaCHIP PG -p NorthSea_67GoodBins -r c -t 12



# if you want to look at HGTs at multiple taxonomic ranks (e.g. phylum, order and genus), 
# you can (and should) skip the gene prediction and blastn steps by providing the "-grouping_only" option.
MetaCHIP PI -i NorthSea_67GoodBins -x fasta -taxon NorthSea_67GoodBins.bac120.summary.tsv -p NorthSea_67GoodBins -r p -t 12 -grouping_only 
MetaCHIP PI -i NorthSea_67GoodBins -x fasta -taxon NorthSea_67GoodBins.bac120.summary.tsv -p NorthSea_67GoodBins -r o -t 12 -grouping_only 
MetaCHIP PI -i NorthSea_67GoodBins -x fasta -taxon NorthSea_67GoodBins.bac120.summary.tsv -p NorthSea_67GoodBins -r g -t 12 -grouping_only 

MetaCHIP BM -p NorthSea_67GoodBins -r p -t 12
MetaCHIP BM -p NorthSea_67GoodBins -r o -t 12
MetaCHIP BM -p NorthSea_67GoodBins -r g -t 12

MetaCHIP PG -p NorthSea_67GoodBins -r p -t 12
MetaCHIP PG -p NorthSea_67GoodBins -r o -t 12
MetaCHIP PG -p NorthSea_67GoodBins -r g -t 12

