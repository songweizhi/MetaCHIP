# working directory
cd /srv/scratch/z5039045/MetaCHIP/BBAY68_69_70/bins


# run Prokka
module load python/3.5.2
module load blast+/2.6.0
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/InPrep1_Prokka.py -genome_folder bins -genome_extension fasta


# run get_homologues
# unload all loaded modules, then
cd /srv/scratch/z5039045/MetaCHIP/BBAY68_69_70/28_good_bins
/srv/scratch/z5039045/Softwares/gh/get_homologues-x86_64-20160712/get_homologues.pl -f 70 -t 3 -S 70 -E 1e-05 -C 70 -G -d input_genomes_gbk


# Test hmmerWrapper
module load python/3.5.2
module load hmmer/3.1b2
module load mafft/7.310
module load fasttree/2.1.7

cd /srv/scratch/z5039045/MetaCHIP/BBAY68_69_70/hmmTree
python3 hmmerWrapper.py 






module load python/3.5.2
module load blast+/2.6.0
cd /srv/scratch/z5039045/MetaCHIP/BBAY68_69_70/21bins
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/InPrep1_Prokka.py -genome_folder Renamed_refined_bins -genome_extension fasta
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/InPrep2_BLAST.py
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/MetaCHIP.py -cfg config.txt



cd /srv/scratch/z5039045/MetaCHIP/m5/6_MetaCHIP_9_million
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/InPrep1_Prokka.py -genome_folder Renamed_refined_bins -genome_extension fasta









############################## Performance on real datasets ##############################

cd /srv/scratch/z5039045/BBAYs/BBAY68_69_70

# combine reads
cat BBAY68_R1_Q20_P.fastq BBAY69_R1_Q20_P.fastq BBAY70_R1_Q20_P.fastq > BBAY68_69_70_R1.fastq
cat BBAY68_R2_Q20_P.fastq BBAY69_R2_Q20_P.fastq BBAY70_R2_Q20_P.fastq > BBAY68_69_70_R2.fastq

# metaSPAdes
module load spades/3.9.0
spades.py --meta -t 12 -1 BBAY68_69_70_R1.fastq -2 BBAY68_69_70_R2.fastq -o BBAY68_69_70_metaSPAdes_k21-127 -k 21,33,55,77,99,127



module load python/2.7.12
module load ete3/3.0.0b36
virtualenv --system-site-packages mypythonenv
. mypythonenv/bin/activate
export PATH=~/anaconda_ete/bin:$PATH;
cd /srv/scratch/z5039045/MetaCHIP/BBAY68_69_70/21bins_renamed
python /srv/scratch/z5039045/Softwares/MetaCHIP/MetaCHIP.py -cfg config.txt
python /srv/scratch/z5039045/Softwares/MetaCHIP/Explicit_Tree.py -cfg config.txt









############################## Performance on 28 good bins ###############################

cd /srv/scratch/z5039045/MetaCHIP/BBAY68_69_70/28_good_bins

# run Prokka
module unload python/2.7.12
module load python/3.5.2
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/InPrep1_Prokka.py -genome_folder 28_good_bins -genome_extension fasta

# run Blast
module unload python/2.7.12
module load python/3.5.2 
python3 /srv/scratch/z5039045/Softwares/MetaCHIP_backup/InPrep2_BLAST.py


# get species tree and grouping input genomes
module unload python/3.5.2 
module load python/2.7.12
module load blast+/2.6.0
module load mafft/7.310
module load hmmer/3.1b2
module load fasttree/2.1.7
python /srv/scratch/z5039045/Softwares/MetaCHIP/Get_species_tree.py -prokka_output prokka_output -hmm /srv/scratch/z5039045/Softwares/MetaCHIP/phylo.hmm 

# run get_homologues
# unload all loaded modules, then
cd /srv/scratch/z5039045/MetaCHIP/BBAY68_69_70/28_good_bins
/srv/scratch/z5039045/Softwares/gh/get_homologues-x86_64-20160712/get_homologues.pl -f 70 -t 3 -S 70 -E 1e-05 -C 70 -G -d input_genomes_gbk


# run MetaCHIP
module unload python/3.5.2 
module load python/2.7.12
module load ete3/3.0.0b36
cd 
virtualenv --system-site-packages mypythonenv
. mypythonenv/bin/activate
export PATH=~/anaconda_ete/bin:$PATH;
cd /srv/scratch/z5039045/MetaCHIP/BBAY68_69_70/28_good_bins
python /srv/scratch/z5039045/Softwares/MetaCHIP/MetaCHIP.py -cfg config.txt
python /srv/scratch/z5039045/Softwares/MetaCHIP/Explicit_Tree.py -cfg config.txt

# COG annotation
module unload python/2.7.12
module load python/3.5.2
module load perl/5.20.1
module load blast+/2.6.0
cd /srv/scratch/z5039045/MetaCHIP/BBAY68_69_70/28_good_bins/output_ip90_al200bp_c70_e1000bp
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in HGT_candidates.fasta -t N
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in HGT_candidates_ET_validated.fasta -t N


# prepare heatplot
python3 ~/PycharmProjects/MetaCHIP/Temp/COG_summary.py -in prokka_output_COG -out COG_annotation000 -in_percent -columns Bin7,Bin96,Bin82,Bin170,Bin71,Bin64,Bin100,Bin129,Bin78,Bin160,Bin99,Bin57,Bin43,Bin74,Bin86,Bin131,Bin162,Bin155,Bin180,Bin25,Bin95,Blast,Tree


# get species tree
module unload python/3.5.2 
module load python/2.7.12


# prepare input

cp -r 82bins_MetaCHIP_wd/prokka_output/


module load python/3.5.2
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/InPrep1_Prokka.py -genome_folder 26bins_2_40 -genome_extension fasta
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/InPrep1_Prokka.py -genome_folder 47bins_5_40 -genome_extension fasta

python3 /srv/scratch/z5039045/Softwares/MetaCHIP/InPrep2
python3 /srv/scratch/z5039045/Softwares/MetaCHIP_backup/InPrep2_BLAST.py


/srv/scratch/z5039045/Softwares/gh/get_homologues-x86_64-20160712/get_homologues.pl -f 70 -t 3 -S 70 -E 1e-05 -C 70 -G -d input_genomes_gbk
/srv/scratch/z5039045/Softwares/gh/get_homologues-x86_64-20160712/get_homologues.pl -f 70 -t 3 -S 70 -E 1e-05 -C 70 -G -d input_genomes_gbk



module unload python/3.5.2 
module load python/2.7.12
module load blast+/2.6.0
module load mafft/7.310
module load hmmer/3.1b2
module load fasttree/2.1.7

python /srv/scratch/z5039045/Softwares/MetaCHIP/Get_species_tree.py -prokka_output prokka_output -hmm /srv/scratch/z5039045/Softwares/MetaCHIP/phylo.hmm 





module unload python/3.5.2 
module load python/2.7.12
module load ete3/3.0.0b36
cd 
virtualenv --system-site-packages mypythonenv
. mypythonenv/bin/activate
export PATH=~/anaconda_ete/bin:$PATH;
cd /srv/scratch/z5039045/MetaCHIP/BBAY68_69_70/47bins_MetaCHIP_wd
python /srv/scratch/z5039045/Softwares/MetaCHIP/Explicit_Tree.py -cfg config.txt


# run test 
module unload python/3.5.2 
module load python/2.7.12
module load ete3/3.0.0b36
cd 
virtualenv --system-site-packages mypythonenv
. mypythonenv/bin/activate
export PATH=~/anaconda_ete/bin:$PATH;
cd /srv/scratch/z5039045/MetaCHIP/BBAY68_69_70/26bins_test
python /srv/scratch/z5039045/Softwares/MetaCHIP/MetaCHIP.py -cfg config.txt
python /srv/scratch/z5039045/Softwares/MetaCHIP/Explicit_Tree.py -cfg config.txt


# COG annotation of tree validated HGT candidates
module load python/3.5.2
module load perl/5.20.1
module load blast+/2.6.0

# 26bin
cd /srv/scratch/z5039045/MetaCHIP/BBAY68_69_70/26bins_MetaCHIP_wd/output_ip90_al200bp_c70_e1000bp
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in HGT_candidates_ET_validated.fasta -t N

# 47bin
cd /srv/scratch/z5039045/MetaCHIP/BBAY68_69_70/47bins_MetaCHIP_wd/output_ip90_al200bp_c70_e1000bp
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in HGT_candidates_ET_validated.fasta -t N


python3 ~/PycharmProjects/MetaCHIP/Temp/COG_summary_boxplot.py -in COG_annotation_26bins -out COG_summary_26bins.csv
python3 ~/PycharmProjects/MetaCHIP/Temp/COG_summary_boxplot.py -in COG_annotation_26bins -out COG_summary_26bins_in_percent.csv -in_percent 
python3 ~/PycharmProjects/MetaCHIP/Temp/COG_summary_boxplot.py -in COG_annotation_47bins -out COG_summary_47bins.csv
python3 ~/PycharmProjects/MetaCHIP/Temp/COG_summary_boxplot.py -in COG_annotation_47bins -out COG_summary_47bins_in_percent.csv -in_percent 




##################################### NorthSea bins ######################################

module unload python/2.7.12
module load python/3.5.2
module load perl/5.20.1
module load blast+/2.6.0

# get quality with universal SCGs
python3 /srv/scratch/z5039045/Scripts/CheckM_runner_universal.py -1 bins_all -pbs -qsub
python3 /srv/scratch/z5039045/Binning_refiner/CheckM_runner.py -1 bins_all -pbs -qsub

# run Prokka
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/InPrep1_Prokka.py -genome_folder bins_37_combined -genome_extension fasta
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/InPrep1_Prokka.py -genome_folder 50bins -genome_extension fasta

# run all_vs_all blastn
python3 /srv/scratch/z5039045/Softwares/MetaCHIP_backup/InPrep2_BLAST.py

# run get homologus
/srv/scratch/z5039045/Softwares/gh/get_homologues-x86_64-20170918/get_homologues.pl -f 70 -t 3 -S 70 -E 1e-05 -C 70 -G -n 6 -d input_genomes_gbk

# get species tree
module unload python/3.5.2 
module load python/2.7.12
module load blast+/2.6.0
module load mafft/7.310
module load hmmer/3.1b2
module load fasttree/2.1.7
python /srv/scratch/z5039045/Softwares/MetaCHIP/Get_species_tree.py -prokka_output prokka_output -hmm /srv/scratch/z5039045/Softwares/MetaCHIP/phylo.hmm 

# get grouping
module load R/3.2.2
R
install.packages('ape')
Rscript /srv/scratch/z5039045/Scripts/genome_grouping.R -t species_tree.newick -c 2 -o grouping_out.txt
Rscript /srv/scratch/z5039045/Scripts/genome_grouping.R -t species_tree.newick -c 5 -o grouping_out.txt
Rscript /srv/scratch/z5039045/Scripts/genome_grouping.R -t species_tree.newick -c 6 -o grouping_out.txt
Rscript /srv/scratch/z5039045/Scripts/genome_grouping.R -t species_tree.newick -c 8 -o grouping_out.txt
Rscript /srv/scratch/z5039045/Scripts/genome_grouping.R -t species_tree.newick -c 9 -o grouping_out.txt
Rscript /srv/scratch/z5039045/Scripts/genome_grouping.R -t species_tree.newick -c 10 -o grouping_out.txt

# turn cluster output into grouping file
python /srv/scratch/z5039045/Scripts/cluster_2_grouping_file.py -c grouping_out.txt -g grouping.txt

# run MetaCHIP
module unload python/3.5.2 
module load python/2.7.12
module load ete3/3.0.0b36
cd 
virtualenv --system-site-packages mypythonenv
. mypythonenv/bin/activate
export PATH=~/anaconda_ete/bin:$PATH;
cd /srv/scratch/z5039045/MetaCHIP/NorthSea_samples/50bins_wd
python /srv/scratch/z5039045/Softwares/MetaCHIP/MetaCHIP.py -g grouping_out.txt -a prokka_output -n all_vs_all_ffn.tab
python /srv/scratch/z5039045/Softwares/MetaCHIP/MetaCHIP.py -g grouping.txt -a prokka_output -n all_vs_all_ffn.tab
python /srv/scratch/z5039045/Softwares/MetaCHIP/MetaCHIP.py -g grouping_out.txt -a prokka_output -n all_vs_all_ffn.tab -l 201

python /srv/scratch/z5039045/Softwares/MetaCHIP/MetaCHIP.py -g grouping_10_modified.txt -a prokka_output -n all_vs_all_ffn.tab


python /srv/scratch/z5039045/Softwares/MetaCHIP/MetaCHIP.py -g grouping_5_manual.txt -a prokka_output -n all_vs_all_ffn.tab -l 201
python /srv/scratch/z5039045/Softwares/MetaCHIP/Explicit_Tree.py -a prokka_output -o orthologs_folder -m /srv/scratch/z5039045/Softwares/MetaCHIP/phylo.hmm 

# COG annotation of predicted HGTs and input genomes
module load python/3.5.2
module load perl/5.20.1
module load blast+/2.6.0
cd /srv/scratch/z5039045/MetaCHIP/NorthSea_samples/47bins_wd/output_ip90_al200bp_c70_e1000bp
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in HGT_candidates.fasta -t N
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in HGT_candidates_ET_validated.fasta -t N
# cd to the output folder of Prokka for each input genomes, than run
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin014.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin019.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin022.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin024.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin025.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin028.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin029.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin034.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin038.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin039.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin040.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin043.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin046.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin050.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin051.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin052.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin053.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin054.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin057.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin062.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin067.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin068.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin069.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin070.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin072.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin074.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin077.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin080.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin082.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin083.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin086.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin089.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin090.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin091.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin092.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin096.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin097.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin099.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin102.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin103.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin104.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin106.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin107.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin111.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin113.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin123.faa -t P &
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in NorthSea_bin125.faa -t P &

# copy COG annotation file to a folder named COG_annotation
cd /srv/scratch/z5039045/MetaCHIP/NorthSea_samples/47bins_wd
mkdir COG_annotation
cp output_ip90_al200bp_c70_e1000bp/HGT_candidates_COG_results/func_stats.txt COG_annotation/Best-match_func_stats.txt &
cp output_ip90_al200bp_c70_e1000bp/HGT_candidates_ET_validated_COG_results/func_stats.txt COG_annotation/Phylogenetic_func_stats.txt &
cp prokka_output/NorthSea_bin014/NorthSea_bin014_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin014_func_stats.txt &
cp prokka_output/NorthSea_bin019/NorthSea_bin019_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin019_func_stats.txt &
cp prokka_output/NorthSea_bin022/NorthSea_bin022_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin022_func_stats.txt &
cp prokka_output/NorthSea_bin024/NorthSea_bin024_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin024_func_stats.txt &
cp prokka_output/NorthSea_bin025/NorthSea_bin025_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin025_func_stats.txt &
cp prokka_output/NorthSea_bin028/NorthSea_bin028_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin028_func_stats.txt &
cp prokka_output/NorthSea_bin029/NorthSea_bin029_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin029_func_stats.txt &
cp prokka_output/NorthSea_bin034/NorthSea_bin034_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin034_func_stats.txt &
cp prokka_output/NorthSea_bin038/NorthSea_bin038_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin038_func_stats.txt &
cp prokka_output/NorthSea_bin039/NorthSea_bin039_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin039_func_stats.txt &
cp prokka_output/NorthSea_bin040/NorthSea_bin040_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin040_func_stats.txt &
cp prokka_output/NorthSea_bin043/NorthSea_bin043_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin043_func_stats.txt &
cp prokka_output/NorthSea_bin046/NorthSea_bin046_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin046_func_stats.txt &
cp prokka_output/NorthSea_bin050/NorthSea_bin050_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin050_func_stats.txt &
cp prokka_output/NorthSea_bin051/NorthSea_bin051_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin051_func_stats.txt &
cp prokka_output/NorthSea_bin052/NorthSea_bin052_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin052_func_stats.txt &
cp prokka_output/NorthSea_bin053/NorthSea_bin053_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin053_func_stats.txt &
cp prokka_output/NorthSea_bin054/NorthSea_bin054_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin054_func_stats.txt &
cp prokka_output/NorthSea_bin057/NorthSea_bin057_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin057_func_stats.txt &
cp prokka_output/NorthSea_bin062/NorthSea_bin062_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin062_func_stats.txt &
cp prokka_output/NorthSea_bin067/NorthSea_bin067_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin067_func_stats.txt &
cp prokka_output/NorthSea_bin068/NorthSea_bin068_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin068_func_stats.txt &
cp prokka_output/NorthSea_bin069/NorthSea_bin069_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin069_func_stats.txt &
cp prokka_output/NorthSea_bin070/NorthSea_bin070_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin070_func_stats.txt &
cp prokka_output/NorthSea_bin072/NorthSea_bin072_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin072_func_stats.txt &
cp prokka_output/NorthSea_bin074/NorthSea_bin074_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin074_func_stats.txt &
cp prokka_output/NorthSea_bin077/NorthSea_bin077_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin077_func_stats.txt &
cp prokka_output/NorthSea_bin080/NorthSea_bin080_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin080_func_stats.txt &
cp prokka_output/NorthSea_bin082/NorthSea_bin082_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin082_func_stats.txt &
cp prokka_output/NorthSea_bin083/NorthSea_bin083_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin083_func_stats.txt &
cp prokka_output/NorthSea_bin086/NorthSea_bin086_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin086_func_stats.txt &
cp prokka_output/NorthSea_bin089/NorthSea_bin089_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin089_func_stats.txt &
cp prokka_output/NorthSea_bin090/NorthSea_bin090_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin090_func_stats.txt &
cp prokka_output/NorthSea_bin091/NorthSea_bin091_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin091_func_stats.txt &
cp prokka_output/NorthSea_bin092/NorthSea_bin092_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin092_func_stats.txt &
cp prokka_output/NorthSea_bin096/NorthSea_bin096_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin096_func_stats.txt &
cp prokka_output/NorthSea_bin097/NorthSea_bin097_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin097_func_stats.txt &
cp prokka_output/NorthSea_bin099/NorthSea_bin099_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin099_func_stats.txt &
cp prokka_output/NorthSea_bin102/NorthSea_bin102_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin102_func_stats.txt &
cp prokka_output/NorthSea_bin103/NorthSea_bin103_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin103_func_stats.txt &
cp prokka_output/NorthSea_bin104/NorthSea_bin104_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin104_func_stats.txt &
cp prokka_output/NorthSea_bin106/NorthSea_bin106_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin106_func_stats.txt &
cp prokka_output/NorthSea_bin107/NorthSea_bin107_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin107_func_stats.txt &
cp prokka_output/NorthSea_bin111/NorthSea_bin111_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin111_func_stats.txt &
cp prokka_output/NorthSea_bin113/NorthSea_bin113_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin113_func_stats.txt &
cp prokka_output/NorthSea_bin123/NorthSea_bin123_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin123_func_stats.txt &
cp prokka_output/NorthSea_bin125/NorthSea_bin125_COG_results/func_stats.txt ./COG_annotation/NorthSea_bin125_func_stats.txt &

# get boxplot for COG annotation (download folder to Desktop)
python3 ~/PycharmProjects/MetaCHIP/Temp/COG_summary_boxplot.py -in COG_annotation -out COG_summary_47bins.csv
python3 ~/PycharmProjects/MetaCHIP/Temp/COG_summary_boxplot.py -in COG_annotation -out COG_summary_47bins_in_percent.csv -in_percent 
python3 ~/PycharmProjects/MetaCHIP/Temp/COG_summary_boxplot.py -in COG_annotation_26bins -out COG_summary_26bins_in_percent.csv -in_percent 

Rscript ~/Rscripts/COG_summary_boxplot.R -i COG_summary_47bins_in_percent.csv -o HGTs_COG



#test
python /srv/scratch/z5039045/MetaCHIP.py -g grouping_out.txt -a prokka_output -n all_vs_all_ffn.tab
python ~/PycharmProjects/MetaCHIP/MetaCHIP.py -g grouping_out.txt -a prokka_output -n all_vs_all_ffn_10000.tab





python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in HGT_candidates.fasta -t P &

