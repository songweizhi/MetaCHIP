
########################## simulate gene transfers with HgtSIM ###########################

python3 HgtSIM.py -t sequences_of_gene_transfers.fasta -i 0 -d distribution_of_transfers.txt -f selected_10_Betaproteobacteria -r 1-0-1-1 -x fna -lf TAGATGAGTGATTAGTTAGTTA -rf TAGATGAGTGATTAGTTAGTTA
python3 HgtSIM.py -t sequences_of_gene_transfers.fasta -i 5 -d distribution_of_transfers.txt -f selected_10_Betaproteobacteria -r 1-0-1-1 -x fna -lf TAGATGAGTGATTAGTTAGTTA -rf TAGATGAGTGATTAGTTAGTTA
python3 HgtSIM.py -t sequences_of_gene_transfers.fasta -i 10 -d distribution_of_transfers.txt -f selected_10_Betaproteobacteria -r 1-0-1-1 -x fna -lf TAGATGAGTGATTAGTTAGTTA -rf TAGATGAGTGATTAGTTAGTTA
python3 HgtSIM.py -t sequences_of_gene_transfers.fasta -i 15 -d distribution_of_transfers.txt -f selected_10_Betaproteobacteria -r 1-0-1-1 -x fna -lf TAGATGAGTGATTAGTTAGTTA -rf TAGATGAGTGATTAGTTAGTTA
python3 HgtSIM.py -t sequences_of_gene_transfers.fasta -i 20 -d distribution_of_transfers.txt -f selected_10_Betaproteobacteria -r 1-0-1-1 -x fna -lf TAGATGAGTGATTAGTTAGTTA -rf TAGATGAGTGATTAGTTAGTTA
python3 HgtSIM.py -t sequences_of_gene_transfers.fasta -i 25 -d distribution_of_transfers.txt -f selected_10_Betaproteobacteria -r 1-0-1-1 -x fna -lf TAGATGAGTGATTAGTTAGTTA -rf TAGATGAGTGATTAGTTAGTTA
python3 HgtSIM.py -t sequences_of_gene_transfers.fasta -i 30 -d distribution_of_transfers.txt -f selected_10_Betaproteobacteria -r 1-0-1-1 -x fna -lf TAGATGAGTGATTAGTTAGTTA -rf TAGATGAGTGATTAGTTAGTTA


#################### Performance assessment without reads simulation #####################

# run Get_clusters.py
python Get_clusters.py -i input_genomes_m0 -x fasta -p M0
python Get_clusters.py -i input_genomes_m5 -x fasta -p M5
python Get_clusters.py -i input_genomes_m10 -x fasta -p M10
python Get_clusters.py -i input_genomes_m15 -x fasta -p M15
python Get_clusters.py -i input_genomes_m20 -x fasta -p M20
python Get_clusters.py -i input_genomes_m25 -x fasta -p M25
python Get_clusters.py -i input_genomes_m30 -x fasta -p M30

# run Best-match.py
python Best-match.py -p M0
python Best-match.py -p M5
python Best-match.py -p M10
python Best-match.py -p M15
python Best-match.py -p M20
python Best-match.py -p M25
python Best-match.py -p M30

# run Phylogenetic.py
python Phylogenetic.py -p M0
python Phylogenetic.py -p M5
python Phylogenetic.py -p M10
python Phylogenetic.py -p M15
python Phylogenetic.py -p M20
python Phylogenetic.py -p M25
python Phylogenetic.py -p M30


##################### Performance assessment with reads simulation #######################

# Step_1_1_GemSIM
module load python/3.5.2
cd /srv/scratch/z5039045/MetaCHIP/m0
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_3_million_m0.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_6_million_m0.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_9_million_m0.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_12_million_m0.txt
cd /srv/scratch/z5039045/MetaCHIP/m5
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_3_million_m5.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_6_million_m5.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_9_million_m5.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_12_million_m5.txt
cd /srv/scratch/z5039045/MetaCHIP/m10
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_3_million_m10.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_6_million_m10.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_9_million_m10.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_12_million_m10.txt
cd /srv/scratch/z5039045/MetaCHIP/m15
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_3_million_m15.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_6_million_m15.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_9_million_m15.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_12_million_m15.txt
cd /srv/scratch/z5039045/MetaCHIP/m20
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_3_million_m20.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_6_million_m20.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_9_million_m20.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_12_million_m20.txt
cd /srv/scratch/z5039045/MetaCHIP/m25
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_3_million_m25.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_6_million_m25.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_9_million_m25.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_12_million_m25.txt
cd /srv/scratch/z5039045/MetaCHIP/m30
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_3_million_m30.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_6_million_m30.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_9_million_m30.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_12_million_m30.txt
cd /srv/scratch/z5039045/MetaCHIP/mNC
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_1_GemSIM.py -cfg config_9_million_mNC.txt

# Step_1_2_QC
module load python/3.5.2
cd /srv/scratch/z5039045/MetaCHIP/m0
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_3_million_m0.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_6_million_m0.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_9_million_m0.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_12_million_m0.txt
cd /srv/scratch/z5039045/MetaCHIP/m5
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_3_million_m5.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_6_million_m5.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_9_million_m5.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_12_million_m5.txt
cd /srv/scratch/z5039045/MetaCHIP/m10
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_3_million_m10.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_6_million_m10.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_9_million_m10.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_12_million_m10.txt
cd /srv/scratch/z5039045/MetaCHIP/m15
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_3_million_m15.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_6_million_m15.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_9_million_m15.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_12_million_m15.txt
cd /srv/scratch/z5039045/MetaCHIP/m20
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_3_million_m20.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_6_million_m20.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_9_million_m20.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_12_million_m20.txt
cd /srv/scratch/z5039045/MetaCHIP/m25
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_3_million_m25.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_6_million_m25.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_9_million_m25.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_12_million_m25.txt
cd /srv/scratch/z5039045/MetaCHIP/m30
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_3_million_m30.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_6_million_m30.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_9_million_m30.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_12_million_m30.txt
cd /srv/scratch/z5039045/MetaCHIP/mNC
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_1_2_QC.py -cfg config_9_million_mNC.txt

# Step_2.1_Assemble_metaSPAdes
module load python/3.5.2
cd /srv/scratch/z5039045/MetaCHIP/m0
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_3_million_m0.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_6_million_m0.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_9_million_m0.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_12_million_m0.txt
cd /srv/scratch/z5039045/MetaCHIP/m5
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_3_million_m5.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_6_million_m5.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_9_million_m5.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_12_million_m5.txt
cd /srv/scratch/z5039045/MetaCHIP/m10
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_3_million_m10.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_6_million_m10.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_9_million_m10.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_12_million_m10.txt
cd /srv/scratch/z5039045/MetaCHIP/m15
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_3_million_m15.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_6_million_m15.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_9_million_m15.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_12_million_m15.txt
cd /srv/scratch/z5039045/MetaCHIP/m20
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_3_million_m20.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_6_million_m20.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_9_million_m20.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_12_million_m20.txt
cd /srv/scratch/z5039045/MetaCHIP/m25
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_3_million_m25.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_6_million_m25.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_9_million_m25.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_12_million_m25.txt
cd /srv/scratch/z5039045/MetaCHIP/m30
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_3_million_m30.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_6_million_m30.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_9_million_m30.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_12_million_m30.txt
cd /srv/scratch/z5039045/MetaCHIP/mNC
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_metaSPAdes.py -cfg config_9_million_mNC.txt

# Step_2.2_Assemble_idba_ud
cd /srv/scratch/z5039045/MetaCHIP/m0
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_3_million_m0.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_6_million_m0.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_9_million_m0.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_12_million_m0.txt
cd /srv/scratch/z5039045/MetaCHIP/m5
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_3_million_m5.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_6_million_m5.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_9_million_m5.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_12_million_m5.txt
cd /srv/scratch/z5039045/MetaCHIP/m10
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_3_million_m10.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_6_million_m10.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_9_million_m10.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_12_million_m10.txt
cd /srv/scratch/z5039045/MetaCHIP/m15
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_3_million_m15.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_6_million_m15.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_9_million_m15.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_12_million_m15.txt
cd /srv/scratch/z5039045/MetaCHIP/m20
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_3_million_m20.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_6_million_m20.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_9_million_m20.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_12_million_m20.txt
cd /srv/scratch/z5039045/MetaCHIP/m25
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_3_million_m25.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_6_million_m25.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_9_million_m25.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_12_million_m25.txt
cd /srv/scratch/z5039045/MetaCHIP/m30
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_3_million_m30.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_6_million_m30.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_9_million_m30.txt
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_12_million_m30.txt
cd /srv/scratch/z5039045/MetaCHIP/mNC
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/assessment_scripts/Step_2_Assemble_idba_ud.py -cfg config_9_million_mNC.txt

# get the number of recovered gene transfers by assemblers from dornor, recipient or both
module load python/3.5.2
module load blast+/2.6.0

python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers_iden100.py -a m0_IDBA_UD_3_million_k20-124.fasta -t input_sequence_mutant_nc_m0.fasta -minf 1000 -d distribution_of_transfers.txt -combined_genomes combined_reference_genomes_m0.fna 
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers_iden100.py -a m0_IDBA_UD_6_million_k20-124.fasta -t input_sequence_mutant_nc_m0.fasta -minf 1000 -d distribution_of_transfers.txt -combined_genomes combined_reference_genomes_m0.fna 
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers_iden100.py -a m0_IDBA_UD_9_million_k20-124.fasta -t input_sequence_mutant_nc_m0.fasta -minf 1000 -d distribution_of_transfers.txt -combined_genomes combined_reference_genomes_m0.fna 
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers_iden100.py -a m0_IDBA_UD_12_million_k20-124.fasta -t input_sequence_mutant_nc_m0.fasta -minf 1000 -d distribution_of_transfers.txt -combined_genomes combined_reference_genomes_m0.fna 
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers_iden100.py -a m0_metaSPAdes_3_million_scaffolds.fasta -t input_sequence_mutant_nc_m0.fasta -minf 1000 -d distribution_of_transfers.txt -combined_genomes combined_reference_genomes_m0.fna 
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers_iden100.py -a m0_metaSPAdes_6_million_scaffolds.fasta -t input_sequence_mutant_nc_m0.fasta -minf 1000 -d distribution_of_transfers.txt -combined_genomes combined_reference_genomes_m0.fna 
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers_iden100.py -a m0_metaSPAdes_9_million_scaffolds.fasta -t input_sequence_mutant_nc_m0.fasta -minf 1000 -d distribution_of_transfers.txt -combined_genomes combined_reference_genomes_m0.fna 
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers_iden100.py -a m0_metaSPAdes_12_million_scaffolds.fasta -t input_sequence_mutant_nc_m0.fasta -minf 1000 -d distribution_of_transfers.txt -combined_genomes combined_reference_genomes_m0.fna 

python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m5_IDBA_UD_3_million_k20-124.fasta -t input_sequence_mutant_nc_m5.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m5_IDBA_UD_6_million_k20-124.fasta -t input_sequence_mutant_nc_m5.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m5_IDBA_UD_9_million_k20-124.fasta -t input_sequence_mutant_nc_m5.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m5_IDBA_UD_12_million_k20-124.fasta -t input_sequence_mutant_nc_m5.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m5_metaSPAdes_3_million_scaffolds.fasta -t input_sequence_mutant_nc_m5.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m5_metaSPAdes_6_million_scaffolds.fasta -t input_sequence_mutant_nc_m5.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m5_metaSPAdes_9_million_scaffolds.fasta -t input_sequence_mutant_nc_m5.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m5_metaSPAdes_12_million_scaffolds.fasta -t input_sequence_mutant_nc_m5.fasta -minf 1000

python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m10_IDBA_UD_3_million_k20-124.fasta -t input_sequence_mutant_nc_m10.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m10_IDBA_UD_6_million_k20-124.fasta -t input_sequence_mutant_nc_m10.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m10_IDBA_UD_9_million_k20-124.fasta -t input_sequence_mutant_nc_m10.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m10_IDBA_UD_12_million_k20-124.fasta -t input_sequence_mutant_nc_m10.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m10_metaSPAdes_3_million_scaffolds.fasta -t input_sequence_mutant_nc_m10.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m10_metaSPAdes_6_million_scaffolds.fasta -t input_sequence_mutant_nc_m10.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m10_metaSPAdes_9_million_scaffolds.fasta -t input_sequence_mutant_nc_m10.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m10_metaSPAdes_12_million_scaffolds.fasta -t input_sequence_mutant_nc_m10.fasta -minf 1000

python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m15_IDBA_UD_3_million_k20-124.fasta -t input_sequence_mutant_nc_m15.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m15_IDBA_UD_6_million_k20-124.fasta -t input_sequence_mutant_nc_m15.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m15_IDBA_UD_9_million_k20-124.fasta -t input_sequence_mutant_nc_m15.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m15_IDBA_UD_12_million_k20-124.fasta -t input_sequence_mutant_nc_m15.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m15_metaSPAdes_3_million_scaffolds.fasta -t input_sequence_mutant_nc_m15.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m15_metaSPAdes_6_million_scaffolds.fasta -t input_sequence_mutant_nc_m15.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m15_metaSPAdes_9_million_scaffolds.fasta -t input_sequence_mutant_nc_m15.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m15_metaSPAdes_12_million_scaffolds.fasta -t input_sequence_mutant_nc_m15.fasta -minf 1000

python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m20_IDBA_UD_3_million_k20-124.fasta -t input_sequence_mutant_nc_m20.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m20_IDBA_UD_6_million_k20-124.fasta -t input_sequence_mutant_nc_m20.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m20_IDBA_UD_9_million_k20-124.fasta -t input_sequence_mutant_nc_m20.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m20_IDBA_UD_12_million_k20-124.fasta -t input_sequence_mutant_nc_m20.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m20_metaSPAdes_3_million_scaffolds.fasta -t input_sequence_mutant_nc_m20.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m20_metaSPAdes_6_million_scaffolds.fasta -t input_sequence_mutant_nc_m20.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m20_metaSPAdes_9_million_scaffolds.fasta -t input_sequence_mutant_nc_m20.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m20_metaSPAdes_12_million_scaffolds.fasta -t input_sequence_mutant_nc_m20.fasta -minf 1000

python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m25_IDBA_UD_3_million_k20-124.fasta -t input_sequence_mutant_nc_m25.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m25_IDBA_UD_6_million_k20-124.fasta -t input_sequence_mutant_nc_m25.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m25_IDBA_UD_9_million_k20-124.fasta -t input_sequence_mutant_nc_m25.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m25_IDBA_UD_12_million_k20-124.fasta -t input_sequence_mutant_nc_m25.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m25_metaSPAdes_3_million_scaffolds.fasta -t input_sequence_mutant_nc_m25.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m25_metaSPAdes_6_million_scaffolds.fasta -t input_sequence_mutant_nc_m25.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m25_metaSPAdes_9_million_scaffolds.fasta -t input_sequence_mutant_nc_m25.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m25_metaSPAdes_12_million_scaffolds.fasta -t input_sequence_mutant_nc_m25.fasta -minf 1000

python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m30_IDBA_UD_3_million_k20-124.fasta -t input_sequence_mutant_nc_m30.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m30_IDBA_UD_6_million_k20-124.fasta -t input_sequence_mutant_nc_m30.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m30_IDBA_UD_9_million_k20-124.fasta -t input_sequence_mutant_nc_m30.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m30_IDBA_UD_12_million_k20-124.fasta -t input_sequence_mutant_nc_m30.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m30_metaSPAdes_3_million_scaffolds.fasta -t input_sequence_mutant_nc_m30.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m30_metaSPAdes_6_million_scaffolds.fasta -t input_sequence_mutant_nc_m30.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m30_metaSPAdes_9_million_scaffolds.fasta -t input_sequence_mutant_nc_m30.fasta -minf 1000
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/scripts/Assembler_recovered_transfers.py -a m30_metaSPAdes_12_million_scaffolds.fasta -t input_sequence_mutant_nc_m30.fasta -minf 1000

python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/get_successfully_binned_transfers_m0.py -a m0_metaSPAdes_9_million_scaffolds.fasta -t input_sequence_mutant_nc_m0.fasta -minf 1000 -d distribution_of_transfers.txt -combined_genomes combined_reference_genomes_m0.fna 
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/get_successfully_binned_transfers.py -a m5_IDBA_UD_9_million_k20-124.fasta -t input_sequence_mutant_nc_m5.fasta -minf 1000 -mutation_level 5
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/get_successfully_binned_transfers.py -a m10_IDBA_UD_9_million_k20-124.fasta -t input_sequence_mutant_nc_m10.fasta -minf 1000 -mutation_level 10
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/get_successfully_binned_transfers.py -a m15_IDBA_UD_9_million_k20-124.fasta -t input_sequence_mutant_nc_m15.fasta -minf 1000 -mutation_level 15
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/get_successfully_binned_transfers.py -a m20_IDBA_UD_9_million_k20-124.fasta -t input_sequence_mutant_nc_m20.fasta -minf 1000 -mutation_level 20
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/get_successfully_binned_transfers.py -a m25_IDBA_UD_9_million_k20-124.fasta -t input_sequence_mutant_nc_m25.fasta -minf 1000 -mutation_level 25
python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/get_successfully_binned_transfers.py -a m30_IDBA_UD_9_million_k20-124.fasta -t input_sequence_mutant_nc_m30.fasta -minf 1000 -mutation_level 30


# rename metaSPAdes contig ID
module load python/3.5.2
cd 2_metaSPAdes_3_million/combined_k21-127
python3 /srv/scratch/z5039045/Scripts/Rename_ctg_id.py -in scaffolds.fasta -out scaffolds_renamed.fasta -splitter _ -keep_number 2
cd 2_metaSPAdes_6_million/combined_k21-127
python3 /srv/scratch/z5039045/Scripts/Rename_ctg_id.py -in scaffolds.fasta -out scaffolds_renamed.fasta -splitter _ -keep_number 2
cd 2_metaSPAdes_9_million/combined_k21-127
python3 /srv/scratch/z5039045/Scripts/Rename_ctg_id.py -in scaffolds.fasta -out scaffolds_renamed.fasta -splitter _ -keep_number 2
cd 2_metaSPAdes_12_million/combined_k21-127
python3 /srv/scratch/z5039045/Scripts/Rename_ctg_id.py -in scaffolds.fasta -out scaffolds_renamed.fasta -splitter _ -keep_number 2

# Step_3_Mapping
module load python/3.5.2
cd /srv/scratch/z5039045/MetaCHIP/m0
python3 /srv/scratch/z5039045/MetaCHIP/m0/Step_3_Mapping_metaSPAdes.py -cfg config_9_million_m0.txt
cd /srv/scratch/z5039045/MetaCHIP/m5
python3 /srv/scratch/z5039045/MetaCHIP/m0/Step_3_Mapping_metaSPAdes.py -cfg config_9_million_m5.txt
cd /srv/scratch/z5039045/MetaCHIP/m10
python3 /srv/scratch/z5039045/MetaCHIP/m0/Step_3_Mapping_metaSPAdes.py -cfg config_9_million_m10.txt
cd /srv/scratch/z5039045/MetaCHIP/m5
python3 /srv/scratch/z5039045/MetaCHIP/Step_3_Mapping_idba_ud.py -cfg config_9_million_m5.txt &
cd /srv/scratch/z5039045/MetaCHIP/m10
python3 /srv/scratch/z5039045/MetaCHIP/Step_3_Mapping_idba_ud.py -cfg config_9_million_m10.txt &
cd /srv/scratch/z5039045/MetaCHIP/m15
python3 /srv/scratch/z5039045/MetaCHIP/Step_3_Mapping_idba_ud.py -cfg config_9_million_m15.txt &
cd /srv/scratch/z5039045/MetaCHIP/m20
python3 /srv/scratch/z5039045/MetaCHIP/Step_3_Mapping_idba_ud.py -cfg config_9_million_m20.txt &
cd /srv/scratch/z5039045/MetaCHIP/m25
python3 /srv/scratch/z5039045/MetaCHIP/Step_3_Mapping_idba_ud.py -cfg config_9_million_m25.txt &
cd /srv/scratch/z5039045/MetaCHIP/m30
python3 /srv/scratch/z5039045/MetaCHIP/Step_3_Mapping_idba_ud.py -cfg config_9_million_m30.txt &
cd /srv/scratch/z5039045/MetaCHIP/mNC
python3 /srv/scratch/z5039045/MetaCHIP/Step_3_Mapping_idba_ud.py -cfg config_9_million_mNC.txt &

# Step_4_1_Binning
module load python/3.5.2
cd /srv/scratch/z5039045/MetaCHIP/m0
python3 /srv/scratch/z5039045/MetaCHIP/m0/Step_4_Binning_metaSPAdes.py -cfg config_9_million_m0.txt
cd /srv/scratch/z5039045/MetaCHIP/m5
python3 /srv/scratch/z5039045/MetaCHIP/m0/Step_4_Binning_metaSPAdes.py -cfg config_9_million_m5.txt
cd /srv/scratch/z5039045/MetaCHIP/m10
python3 /srv/scratch/z5039045/MetaCHIP/m0/Step_4_Binning_metaSPAdes.py -cfg config_9_million_m10.txt
cd /srv/scratch/z5039045/MetaCHIP/m5
python3 /srv/scratch/z5039045/MetaCHIP/Step_4_Binning_idba_ud.py -cfg config_9_million_m5.txt
cd /srv/scratch/z5039045/MetaCHIP/m10
python3 /srv/scratch/z5039045/MetaCHIP/Step_4_Binning_idba_ud.py -cfg config_9_million_m10.txt
cd /srv/scratch/z5039045/MetaCHIP/m15
python3 /srv/scratch/z5039045/MetaCHIP/Step_4_Binning_idba_ud.py -cfg config_9_million_m15.txt
cd /srv/scratch/z5039045/MetaCHIP/m20
python3 /srv/scratch/z5039045/MetaCHIP/Step_4_Binning_idba_ud.py -cfg config_9_million_m20.txt
cd /srv/scratch/z5039045/MetaCHIP/m25
python3 /srv/scratch/z5039045/MetaCHIP/Step_4_Binning_idba_ud.py -cfg config_9_million_m25.txt
cd /srv/scratch/z5039045/MetaCHIP/m30
python3 /srv/scratch/z5039045/MetaCHIP/Step_4_Binning_idba_ud.py -cfg config_9_million_m30.txt
cd /srv/scratch/z5039045/MetaCHIP/mNC
python3 /srv/scratch/z5039045/MetaCHIP/Step_4_Binning_idba_ud.py -cfg config_9_million_mNC.txt

# Step_4_2_Binning_refiner
module load python/3.5.2
module load R/3.2.2

cd /srv/scratch/z5039045/MetaCHIP/m0/4_Binning_9_million
mv *_*_*mer_*_cov MyCC
mv *_MetaBAT MetaBAT
python3 /srv/scratch/z5039045/Binning_refiner/Binning_refiner.py -1 MetaBAT -2 MyCC 

cd /srv/scratch/z5039045/MetaCHIP/m5/4_Binning_9_million
mv *_*_*mer_*_cov MyCC
mv *_MetaBAT MetaBAT
python3 /srv/scratch/z5039045/Binning_refiner/Binning_refiner.py -1 MetaBAT -2 MyCC 

cd /srv/scratch/z5039045/MetaCHIP/m10/4_Binning_9_million
mv *_*_*mer_*_cov MyCC
mv *_MetaBAT MetaBAT
python3 /srv/scratch/z5039045/Binning_refiner/Binning_refiner.py -1 MetaBAT -2 MyCC 

cd /srv/scratch/z5039045/MetaCHIP/m15/4_Binning_9_million
mv *_*_*mer_*_cov MyCC
mv *_MetaBAT MetaBAT
python3 /srv/scratch/z5039045/Binning_refiner/Binning_refiner.py -1 MetaBAT -2 MyCC 

cd /srv/scratch/z5039045/MetaCHIP/m20/4_Binning_9_million
mv *_*_*mer_*_cov MyCC
mv *_MetaBAT MetaBAT
python3 /srv/scratch/z5039045/Binning_refiner/Binning_refiner.py -1 MetaBAT -2 MyCC 

cd /srv/scratch/z5039045/MetaCHIP/m25/4_Binning_9_million
mv *_*_*mer_*_cov MyCC
mv *_MetaBAT MetaBAT
python3 /srv/scratch/z5039045/Binning_refiner/Binning_refiner.py -1 MetaBAT -2 MyCC 

cd /srv/scratch/z5039045/MetaCHIP/m30/4_Binning_9_million
mv *_*_*mer_*_cov MyCC
mv *_MetaBAT MetaBAT
python3 /srv/scratch/z5039045/Binning_refiner/Binning_refiner.py -1 MetaBAT -2 MyCC 

cd /srv/scratch/z5039045/MetaCHIP/mNC/4_Binning_9_million
mv *_*_*mer_*_cov MyCC
mv *_MetaBAT MetaBAT
python3 /srv/scratch/z5039045/Binning_refiner/Binning_refiner.py -1 MetaBAT -2 MyCC 


# Step_5_get_correlation
module load python/3.5.2
module load blast+/2.6.0
module load R/3.2.2

# m0
cd /srv/scratch/z5039045/MetaCHIP/m0
mkdir 5_get_correlations_9_million
cp -r 4_Binning_9_million/outputs/Refined 5_get_correlations_9_million
cp -r input_genomes_m0 5_get_correlations_9_million/REF
cd 5_get_correlations_9_million
python3 /srv/scratch/z5039045/Scripts/get_correlations.py -1 Refined -2 REF

cd /srv/scratch/z5039045/MetaCHIP/m5
mkdir 5_get_correlations_9_million
cp -r 4_Binning_9_million/outputs/Refined 5_get_correlations_9_million
cp -r input_genomes_m5 5_get_correlations_9_million/REF
cd 5_get_correlations_9_million
python3 /srv/scratch/z5039045/Scripts/get_correlations.py -1 Refined -2 REF

cd /srv/scratch/z5039045/MetaCHIP/m10
mkdir 5_get_correlations_9_million
cp -r 4_Binning_9_million/outputs/Refined 5_get_correlations_9_million
cp -r input_genomes_m10 5_get_correlations_9_million/REF
cd 5_get_correlations_9_million
python3 /srv/scratch/z5039045/Scripts/get_correlations.py -1 Refined -2 REF

cd /srv/scratch/z5039045/MetaCHIP/m15
mkdir 5_get_correlations_9_million
cp -r 4_Binning_9_million/outputs/Refined 5_get_correlations_9_million
cp -r input_genomes_m15 5_get_correlations_9_million/REF
cd 5_get_correlations_9_million
python3 /srv/scratch/z5039045/Scripts/get_correlations.py -1 Refined -2 REF

cd /srv/scratch/z5039045/MetaCHIP/m20
mkdir 5_get_correlations_9_million
cp -r 4_Binning_9_million/outputs/Refined 5_get_correlations_9_million
cp -r input_genomes_m20 5_get_correlations_9_million/REF
cd 5_get_correlations_9_million
python3 /srv/scratch/z5039045/Scripts/get_correlations.py -1 Refined -2 REF

cd /srv/scratch/z5039045/MetaCHIP/m25
mkdir 5_get_correlations_9_million
cp -r 4_Binning_9_million/outputs/Refined 5_get_correlations_9_million
cp -r input_genomes_m25 5_get_correlations_9_million/REF
cd 5_get_correlations_9_million
python3 /srv/scratch/z5039045/Scripts/get_correlations.py -1 Refined -2 REF

cd /srv/scratch/z5039045/MetaCHIP/m30
mkdir 5_get_correlations_9_million
cp -r 4_Binning_9_million/outputs/Refined 5_get_correlations_9_million
cp -r input_genomes_m30 5_get_correlations_9_million/REF
cd 5_get_correlations_9_million
python3 /srv/scratch/z5039045/Scripts/get_correlations.py -1 Refined -2 REF

# renamed refined bins
cd /srv/scratch/z5039045/MetaCHIP/m0/5_get_correlations_9_million
cp -r Refined Renamed_refined_bins
cd Renamed_refined_bins
cat Refined_24.fasta > AAM.fasta
cat Refined_11.fasta > AKV.fasta
cat Refined_4.fasta > AMAC.fasta
cat Refined_17.fasta > AMAU.fasta
cat Refined_18.fasta > AMS.fasta
cat Refined_20.fasta > ARL.fasta
cat Refined_21.fasta > ARS.fasta
cat Refined_7.fasta Refined_8.fasta > ASJ.fasta
cat Refined_22.fasta > ASN.fasta
cat Refined_5.fasta > ATM.fasta
cat Refined_2.fasta Refined_3.fasta > BAD.fasta
cat Refined_6.fasta > BDS.fasta
cat Refined_10.fasta > BGC.fasta
cat Refined_15.fasta > BHS.fasta
cat Refined_12.fasta > BNM.fasta
cat Refined_1.fasta Refined_16.fasta > BRT.fasta
cat Refined_14.fasta > BSA.fasta
cat Refined_9.fasta > BSD.fasta
cat Refined_13.fasta > BSL.fasta
cat Refined_23.fasta > BTK.fasta

# m5 get bin correlations
cd /srv/scratch/z5039045/MetaCHIP/m5
mkdir 5_get_correlations_9_million
cp -r 4_Binning_9_million/outputs/Refined 5_get_correlations_9_million
cp -r input_genomes_m5 5_get_correlations_9_million/REF
cd 5_get_correlations_9_million
python3 /srv/scratch/z5039045/Scripts/get_correlations.py -1 Refined -2 REF

# m10 get bin correlations
cd /srv/scratch/z5039045/MetaCHIP/m10
mkdir 5_get_correlations_9_million
cp -r 4_Binning_9_million/outputs/Refined 5_get_correlations_9_million
cp -r input_genomes_m10 5_get_correlations_9_million/REF
cd 5_get_correlations_9_million
python3 /srv/scratch/z5039045/Scripts/get_correlations.py -1 Refined -2 REF

# renamed refined bins
cd /srv/scratch/z5039045/MetaCHIP/m5/5_get_correlations_9_million
cp -r Refined Renamed_refined_bins
cd Renamed_refined_bins
cat Refined_16.fasta Refined_5.fasta > AAM.fasta &
cat Refined_10.fasta > AKV.fasta &
cat Refined_24.fasta > AMAC.fasta &
cat Refined_13.fasta Refined_14.fasta > AMAU.fasta &
cat Refined_18.fasta > AMS.fasta &
cat Refined_20.fasta > ARL.fasta &
cat Refined_22.fasta > ARS.fasta &
cat Refined_6.fasta > ASJ.fasta &
cat Refined_23.fasta > ASN.fasta &
cat Refined_3.fasta > ATM.fasta &
cat Refined_15.fasta Refined_2.fasta > BAD.fasta &
cat Refined_4.fasta > BDS.fasta &
cat Refined_9.fasta > BGC.fasta &
cat Refined_17.fasta > BHS.fasta &
cat Refined_7.fasta > BNM.fasta &
cat Refined_1.fasta > BRT.fasta &
cat Refined_11.fasta > BSA.fasta &
cat Refined_8.fasta > BSD.fasta &
cat Refined_12.fasta > BSL.fasta &
cat Refined_21.fasta > BTK.fasta &
rm Refined_*.fasta

cd /srv/scratch/z5039045/MetaCHIP/m10/5_get_correlations_9_million
cp -r Refined Renamed_refined_bins
cd Renamed_refined_bins
cat Refined_5.fasta > AAM.fasta &
cat Refined_11.fasta > AKV.fasta &
cat Refined_3.fasta > AMAC.fasta &
cat Refined_13.fasta > AMAU.fasta &
cat Refined_18.fasta > AMS.fasta &
cat Refined_16.fasta > ARL.fasta &
cat Refined_19.fasta > ARS.fasta &
cat Refined_6.fasta > ASJ.fasta &
cat Refined_15.fasta > ASN.fasta &
cat Refined_20.fasta > ATM.fasta &
cat Refined_1.fasta > BAD.fasta &
cat Refined_4.fasta > BDS.fasta &
cat Refined_10.fasta > BGC.fasta &
cat Refined_14.fasta > BHS.fasta &
cat Refined_8.fasta > BNM.fasta &
cat Refined_2.fasta > BRT.fasta &
cat Refined_12.fasta > BSA.fasta &
cat Refined_7.fasta > BSD.fasta &
cat Refined_9.fasta > BSL.fasta &
cat Refined_17.fasta > BTK.fasta &
rm Refined_*.fasta

cd /srv/scratch/z5039045/MetaCHIP/m15/5_get_correlations_9_million
cp -r Refined Renamed_refined_bins
cd Renamed_refined_bins
cat Refined_6.fasta > AAM.fasta &
cat Refined_11.fasta > AKV.fasta &
cat Refined_3.fasta > AMAC.fasta &
cat Refined_1.fasta > AMAU.fasta &
cat Refined_18.fasta > AMS.fasta &
cat Refined_16.fasta > ARL.fasta &
cat Refined_19.fasta > ARS.fasta &
cat Refined_5.fasta > ASJ.fasta &
cat Refined_14.fasta > ASN.fasta &
cat Refined_2.fasta > ATM.fasta &
cat Refined_15.fasta > BAD.fasta &
cat Refined_4.fasta > BDS.fasta &
cat Refined_10.fasta > BGC.fasta &
cat Refined_12.fasta > BHS.fasta &
cat Refined_8.fasta > BNM.fasta &
cat Refined_20.fasta > BRT.fasta &
cat Refined_13.fasta > BSA.fasta &
cat Refined_7.fasta > BSD.fasta &
cat Refined_9.fasta > BSL.fasta &
cat Refined_17.fasta > BTK.fasta &
rm Refined_*.fasta

cd /srv/scratch/z5039045/MetaCHIP/m20/5_get_correlations_9_million
cp -r Refined Renamed_refined_bins
cd Renamed_refined_bins
cat Refined_5.fasta > AAM.fasta &
cat Refined_11.fasta > AKV.fasta &
cat Refined_3.fasta > AMAC.fasta &
cat Refined_1.fasta > AMAU.fasta &
cat Refined_19.fasta > AMS.fasta &
cat Refined_17.fasta > ARL.fasta &
cat Refined_20.fasta > ARS.fasta &
cat Refined_6.fasta > ASJ.fasta &
cat Refined_16.fasta > ASN.fasta &
cat Refined_2.fasta > ATM.fasta &
cat Refined_15.fasta > BAD.fasta &
cat Refined_4.fasta > BDS.fasta &
cat Refined_10.fasta > BGC.fasta &
cat Refined_12.fasta Refined_13.fasta > BHS.fasta &
cat Refined_9.fasta > BNM.fasta &
cat Refined_21.fasta > BRT.fasta &
cat Refined_14.fasta > BSA.fasta &
cat Refined_8.fasta > BSD.fasta &
cat Refined_7.fasta > BSL.fasta &
cat Refined_18.fasta > BTK.fasta &
rm Refined_*.fasta

cd /srv/scratch/z5039045/MetaCHIP/m25/5_get_correlations_9_million
cp -r Refined Renamed_refined_bins
cd Renamed_refined_bins
cat Refined_5.fasta > AAM.fasta &
cat Refined_11.fasta > AKV.fasta &
cat Refined_3.fasta > AMAC.fasta &
cat Refined_1.fasta > AMAU.fasta &
cat Refined_19.fasta > AMS.fasta &
cat Refined_16.fasta > ARL.fasta &
cat Refined_17.fasta > ARS.fasta &
cat Refined_6.fasta > ASJ.fasta &
cat Refined_14.fasta > ASN.fasta &
cat Refined_2.fasta > ATM.fasta &
cat Refined_15.fasta > BAD.fasta &
cat Refined_4.fasta > BDS.fasta &
cat Refined_10.fasta > BGC.fasta &
cat Refined_12.fasta > BHS.fasta &
cat Refined_9.fasta > BNM.fasta &
cat Refined_20.fasta > BRT.fasta &
cat Refined_13.fasta > BSA.fasta &
cat Refined_7.fasta > BSD.fasta &
cat Refined_8.fasta > BSL.fasta &
cat Refined_18.fasta > BTK.fasta &
rm Refined_*.fasta

cd /srv/scratch/z5039045/MetaCHIP/m30/5_get_correlations_9_million
cp -r Refined Renamed_refined_bins
cd Renamed_refined_bins
cat Refined_2.fasta > AAM.fasta &
cat Refined_11.fasta > AKV.fasta &
cat Refined_4.fasta > AMAC.fasta &
cat Refined_13.fasta > AMAU.fasta &
cat Refined_16.fasta > AMS.fasta &
cat Refined_17.fasta > ARL.fasta &
cat Refined_20.fasta > ARS.fasta &
cat Refined_6.fasta > ASJ.fasta &
cat Refined_1.fasta > ASN.fasta &
cat Refined_3.fasta > ATM.fasta &
cat Refined_19.fasta > BAD.fasta &
cat Refined_5.fasta > BDS.fasta &
cat Refined_10.fasta > BGC.fasta &
cat Refined_15.fasta > BHS.fasta &
cat Refined_8.fasta > BNM.fasta &
cat Refined_21.fasta > BRT.fasta &
cat Refined_12.fasta > BSA.fasta &
cat Refined_7.fasta > BSD.fasta &
cat Refined_9.fasta > BSL.fasta &
cat Refined_18.fasta > BTK.fasta &
rm Refined_*.fasta


# Get the precision and recall of refined genome bins with Evaluate.py 
module unload intel/11.1.080
module unload python/3.5.2
module load intel/16.0.1.150
module load perl/5.20.1
module load python/2.7.10
module load cd-hit/4.6.4
module load prodigal/2.6.3
module load parallel/20160222
module load hmmer/3.1b2
module load barrnap/0.7
module load mycc/20150710
module load blast+/2.6.0

cat m0/input_genomes_m0/*.fna > m0/3_Mapping_9_million/combined_input_genomes.fna
cd /srv/scratch/z5039045/MetaCHIP/m0/3_Mapping_9_million
blastn -query scaffolds_k21-127_lt2000.fa -subject combined_input_genomes.fna -out m0_ctg_assignment.tab -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" &
cd /srv/scratch/z5039045/MetaCHIP/m0/4_Binning_9_million/outputs/Refined
Evaluate.py /srv/scratch/z5039045/MetaCHIP/qsub_blast_precision_recall/m0_ctg_assignment_id_uniq.tab /srv/scratch/z5039045/MetaCHIP/m0/3_Mapping_9_million/scaffolds_k21-127_lt2000.fa fasta 
No. of reference genomes: 20
No. of bins in evaluation: 24
No. of sequences assigned reference: 1049
No. of binned sequences: 1002
Precision: 0.989643, 0.997291
Sensitivity: 0.735939, 0.894861

cat m5/input_genomes_m5/*.fna > m5/3_Mapping_9_million/combined_input_genomes.fna
cd /srv/scratch/z5039045/MetaCHIP/m5/3_Mapping_9_million
blastn -query scaffolds_k21-127_lt2000.fa -subject combined_input_genomes.fna -out m5_ctg_assignment.tab -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" &
cd /srv/scratch/z5039045/MetaCHIP/m5/4_Binning_9_million/outputs/Refined
Evaluate.py /srv/scratch/z5039045/MetaCHIP/qsub_blast_precision_recall/m5_ctg_assignment_id_uniq.tab /srv/scratch/z5039045/MetaCHIP/m5/3_Mapping_9_million/scaffold_k20-124_lt2500.fa fasta 
No. of reference genomes: 20
No. of bins in evaluation: 23
No. of sequences assigned reference: 676
No. of binned sequences: 620
Precision: 0.991039, 0.999612
Sensitivity: 0.754438, 0.931990

cat m10/input_genomes_m10/*.fna > m10/3_Mapping_9_million/combined_input_genomes.fna
cd /srv/scratch/z5039045/MetaCHIP/m10/3_Mapping_9_million
blastn -query scaffolds_k21-127_lt2000.fa -subject combined_input_genomes.fna -out m10_ctg_assignment.tab -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" &
cd /srv/scratch/z5039045/MetaCHIP/m10/4_Binning_9_million/outputs/Refined
Evaluate.py /srv/scratch/z5039045/MetaCHIP/qsub_blast_precision_recall/m10_ctg_assignment_id_uniq.tab /srv/scratch/z5039045/MetaCHIP/m10/3_Mapping_9_million/scaffold_k20-124_lt2500.fa fasta 
No. of reference genomes: 20
No. of bins in evaluation: 20
No. of sequences assigned reference: 561
No. of binned sequences: 481
Precision: 0.988372, 0.999467
Sensitivity: 0.757576, 0.959221

cat m15/input_genomes_m15/*.fna > m15/3_Mapping_9_million/combined_input_genomes.fna
cd /srv/scratch/z5039045/MetaCHIP/m15/3_Mapping_9_million
blastn -query scaffolds_k21-127_lt2000.fa -subject combined_input_genomes.fna -out m15_ctg_assignment.tab -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" &
cd /srv/scratch/z5039045/MetaCHIP/m15/4_Binning_9_million/outputs/Refined
Evaluate.py /srv/scratch/z5039045/MetaCHIP/qsub_blast_precision_recall/m15_ctg_assignment_id_uniq.tab /srv/scratch/z5039045/MetaCHIP/m15/3_Mapping_9_million/scaffold_k20-124_lt2500.fa fasta 
No. of reference genomes: 20
No. of bins in evaluation: 20
No. of sequences assigned reference: 561
No. of binned sequences: 491
Precision: 0.995536, 0.999669
Sensitivity: 0.795009, 0.964600

cat m20/input_genomes_m20/*.fna > m20/3_Mapping_9_million/combined_input_genomes.fna
cd /srv/scratch/z5039045/MetaCHIP/m20/3_Mapping_9_million
blastn -query scaffolds_k21-127_lt2000.fa -subject combined_input_genomes.fna -out m20_ctg_assignment.tab -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" &
cd /srv/scratch/z5039045/MetaCHIP/m20/4_Binning_9_million/outputs/Refined
Evaluate.py /srv/scratch/z5039045/MetaCHIP/qsub_blast_precision_recall/m20_ctg_assignment_id_uniq.tab /srv/scratch/z5039045/MetaCHIP/m20/3_Mapping_9_million/scaffold_k20-124_lt2500.fa fasta 
No. of reference genomes: 20
No. of bins in evaluation: 21
No. of sequences assigned reference: 544
No. of binned sequences: 483
Precision: 0.988345, 0.999291
Sensitivity: 0.768382, 0.954067

cat m25/input_genomes_m25/*.fna > m25/3_Mapping_9_million/combined_input_genomes.fna
cd /srv/scratch/z5039045/MetaCHIP/m25/3_Mapping_9_million
blastn -query scaffolds_k21-127_lt2000.fa -subject combined_input_genomes.fna -out m25_ctg_assignment.tab -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" &
cd /srv/scratch/z5039045/MetaCHIP/m25/4_Binning_9_million/outputs/Refined
Evaluate.py /srv/scratch/z5039045/MetaCHIP/qsub_blast_precision_recall/m25_ctg_assignment_id_uniq.tab /srv/scratch/z5039045/MetaCHIP/m25/3_Mapping_9_million/scaffold_k20-124_lt2500.fa fasta 
No. of reference genomes: 20
No. of bins in evaluation: 20
No. of sequences assigned reference: 537
No. of binned sequences: 472
Precision: 0.995227, 0.999675
Sensitivity: 0.776536, 0.964525

cat m30/input_genomes_m30/*.fna > m30/3_Mapping_9_million/combined_input_genomes.fna
cd /srv/scratch/z5039045/MetaCHIP/m30/3_Mapping_9_million
blastn -query scaffolds_k21-127_lt2000.fa -subject combined_input_genomes.fna -out m30_ctg_assignment.tab -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" &
cd /srv/scratch/z5039045/MetaCHIP/m30/4_Binning_9_million/outputs/Refined
Evaluate.py /srv/scratch/z5039045/MetaCHIP/qsub_blast_precision_recall/m30_ctg_assignment_id_uniq.tab /srv/scratch/z5039045/MetaCHIP/m30/3_Mapping_9_million/scaffold_k20-124_lt2500.fa fasta 
No. of reference genomes: 20
No. of bins in evaluation: 21
No. of sequences assigned reference: 559
No. of binned sequences: 495
Precision: 0.997763, 0.999959
Sensitivity: 0.774597, 0.963484


# Step_6_run_MetaCHIP

mkdir 6_MetaCHIP_9_million
mv 5_get_correlations_9_million/Renamed_refined_bins 6_MetaCHIP_9_million/

cd /srv/scratch/z5039045/MetaCHIP/m5
mkdir 6_MetaCHIP_9_million
mv 5_get_correlations_9_million/Renamed_refined_bins 6_MetaCHIP_9_million/

cd /srv/scratch/z5039045/MetaCHIP/m10
mkdir 6_MetaCHIP_9_million
mv 5_get_correlations_9_million/Renamed_refined_bins 6_MetaCHIP_9_million/

cd /srv/scratch/z5039045/MetaCHIP/m15
mkdir 6_MetaCHIP_9_million
mv 5_get_correlations_9_million/Renamed_refined_bins 6_MetaCHIP_9_million/

cd /srv/scratch/z5039045/MetaCHIP/m20
mkdir 6_MetaCHIP_9_million
mv 5_get_correlations_9_million/Renamed_refined_bins 6_MetaCHIP_9_million/

cd /srv/scratch/z5039045/MetaCHIP/m25
mkdir 6_MetaCHIP_9_million
mv 5_get_correlations_9_million/Renamed_refined_bins 6_MetaCHIP_9_million/

cd /srv/scratch/z5039045/MetaCHIP/m30
mkdir 6_MetaCHIP_9_million
mv 5_get_correlations_9_million/Renamed_refined_bins 6_MetaCHIP_9_million/


# Step_6 prepare input files and run MetaCHIP

# Get_clusters.py
python Get_clusters.py -i input_genomes_m0 -x fasta -p M0
python Get_clusters.py -i input_genomes_m5 -x fasta -p M5
python Get_clusters.py -i input_genomes_m10 -x fasta -p M10
python Get_clusters.py -i input_genomes_m15 -x fasta -p M15
python Get_clusters.py -i input_genomes_m20 -x fasta -p M20
python Get_clusters.py -i input_genomes_m25 -x fasta -p M25
python Get_clusters.py -i input_genomes_m30 -x fasta -p M30

# Best-match.py
python Best-match.py -p M0
python Best-match.py -p M5
python Best-match.py -p M10
python Best-match.py -p M15
python Best-match.py -p M20
python Best-match.py -p M25
python Best-match.py -p M30

# Phylogenetic.py
python Phylogenetic.py -p M0
python Phylogenetic.py -p M5
python Phylogenetic.py -p M10
python Phylogenetic.py -p M15
python Phylogenetic.py -p M20
python Phylogenetic.py -p M25
python Phylogenetic.py -p M30


# Step_7_get_gene_correlations
module load python/3.5.2
module load blast+/2.6.0
module load R/3.2.2

cd /srv/scratch/z5039045/MetaCHIP/m0
mkdir 7_get_gene_correlations_9_million
cp 6_MetaCHIP_9_million/output_ip90_al200bp_c70_e1000bp/HGT_candidates_with_direction.txt 7_get_gene_correlations_9_million
cp 6_MetaCHIP_9_million/output_ip90_al200bp_c70_e1000bp/HGT_candidates_ET_validated.txt 7_get_gene_correlations_9_million
cp 6_MetaCHIP_9_million/combined.ffn 7_get_gene_correlations_9_million
cd /srv/scratch/z5039045/MetaCHIP/m0/7_get_gene_correlations_9_million
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations_iden100.py -predicted_HGTs HGT_candidates_with_direction.txt -t input_sequence_mutant_nc_m0.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn 
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations_iden100.py -predicted_HGTs HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m0.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn 
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations_iden100.py -predicted_HGTs m0_e1000_HGT_candidates.txt -t input_sequence_mutant_nc_m0.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn 
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations_iden100.py -predicted_HGTs m0_e2000_HGT_candidates.txt -t input_sequence_mutant_nc_m0.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn 
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations_iden100.py -predicted_HGTs m0_e5000_HGT_candidates.txt -t input_sequence_mutant_nc_m0.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn 
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations_iden100.py -predicted_HGTs m0_e10000_HGT_candidates.txt -t input_sequence_mutant_nc_m0.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn 
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations_iden100.py -predicted_HGTs m0_e50000_HGT_candidates.txt -t input_sequence_mutant_nc_m0.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn 
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations_iden100.py -predicted_HGTs m0_e1000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m0.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn 
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations_iden100.py -predicted_HGTs m0_e2000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m0.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn 
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations_iden100.py -predicted_HGTs m0_e5000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m0.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn 
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations_iden100.py -predicted_HGTs m0_e10000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m0.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn 
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations_iden100.py -predicted_HGTs m0_e50000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m0.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn 

cd /srv/scratch/z5039045/MetaCHIP/m5
mkdir 7_get_gene_correlations_9_million
cp 6_MetaCHIP_9_million/output_ip90%_al200bp_c70%_e500bp/HGT_candidates_with_direction.txt 7_get_gene_correlations_9_million
cp 6_MetaCHIP_9_million/output_ip90_al200bp_c70_e1000bp/HGT_candidates_ET_validated.txt 7_get_gene_correlations_9_million
cp 6_MetaCHIP_9_million/combined.ffn 7_get_gene_correlations_9_million
cp ../assemblies/input_sequence_mutant_nc_m5.fasta 7_get_gene_correlations_9_million/
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs HGT_candidates_with_direction.txt -t input_sequence_mutant_nc_m5.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m5.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m5_e1000_HGT_candidates.txt -t input_sequence_mutant_nc_m5.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m5_e2000_HGT_candidates.txt -t input_sequence_mutant_nc_m5.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m5_e5000_HGT_candidates.txt -t input_sequence_mutant_nc_m5.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m5_e10000_HGT_candidates.txt -t input_sequence_mutant_nc_m5.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m5_e50000_HGT_candidates.txt -t input_sequence_mutant_nc_m5.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m5_e1000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m5.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m5_e2000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m5.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m5_e5000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m5.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m5_e10000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m5.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m5_e50000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m5.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn

cd /srv/scratch/z5039045/MetaCHIP/m10
mkdir 7_get_gene_correlations_9_million
cp 6_MetaCHIP_9_million/output_ip90%_al200bp_c70%_e500bp/HGT_candidates_with_direction.txt 7_get_gene_correlations_9_million
cp 6_MetaCHIP_9_million/output_ip90_al200bp_c70_e1000bp/HGT_candidates_ET_validated.txt 7_get_gene_correlations_9_million
cp 6_MetaCHIP_9_million/combined.ffn 7_get_gene_correlations_9_million
cp ../assemblies/input_sequence_mutant_nc_m10.fasta 7_get_gene_correlations_9_million/
cd 7_get_gene_correlations_9_million
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs HGT_candidates_with_direction.txt -t input_sequence_mutant_nc_m10.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m10.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m10_e1000_HGT_candidates.txt -t input_sequence_mutant_nc_m10.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m10_e2000_HGT_candidates.txt -t input_sequence_mutant_nc_m10.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m10_e5000_HGT_candidates.txt -t input_sequence_mutant_nc_m10.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m10_e10000_HGT_candidates.txt -t input_sequence_mutant_nc_m10.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m10_e50000_HGT_candidates.txt -t input_sequence_mutant_nc_m10.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m10_e1000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m10.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m10_e2000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m10.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m10_e5000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m10.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m10_e10000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m10.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m10_e50000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m10.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn

cd /srv/scratch/z5039045/MetaCHIP/m15
mkdir 7_get_gene_correlations_9_million
cp 6_MetaCHIP_9_million/output_ip90%_al200bp_c70%_e500bp/HGT_candidates_with_direction.txt 7_get_gene_correlations_9_million
cp 6_MetaCHIP_9_million/output_ip90_al200bp_c70_e1000bp/HGT_candidates_ET_validated.txt 7_get_gene_correlations_9_million
cp 6_MetaCHIP_9_million/combined.ffn 7_get_gene_correlations_9_million
cp ../assemblies/input_sequence_mutant_nc_m15.fasta 7_get_gene_correlations_9_million/
cd 7_get_gene_correlations_9_million
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs HGT_candidates_with_direction.txt -t input_sequence_mutant_nc_m15.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m15.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m15_e1000_HGT_candidates.txt -t input_sequence_mutant_nc_m15.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m15_e2000_HGT_candidates.txt -t input_sequence_mutant_nc_m15.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m15_e5000_HGT_candidates.txt -t input_sequence_mutant_nc_m15.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m15_e10000_HGT_candidates.txt -t input_sequence_mutant_nc_m15.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m15_e50000_HGT_candidates.txt -t input_sequence_mutant_nc_m15.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m15_e1000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m15.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m15_e2000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m15.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m15_e5000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m15.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m15_e10000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m15.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m15_e50000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m15.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn


cd /srv/scratch/z5039045/MetaCHIP/m20
mkdir 7_get_gene_correlations_9_million
cp 6_MetaCHIP_9_million/output_ip90%_al200bp_c70%_e500bp/HGT_candidates_with_direction.txt 7_get_gene_correlations_9_million
cp 6_MetaCHIP_9_million/output_ip90_al200bp_c70_e1000bp/HGT_candidates_ET_validated.txt 7_get_gene_correlations_9_million
cp 6_MetaCHIP_9_million/combined.ffn 7_get_gene_correlations_9_million
cp ../assemblies/input_sequence_mutant_nc_m20.fasta 7_get_gene_correlations_9_million/
cd 7_get_gene_correlations_9_million
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs HGT_candidates_with_direction.txt -t input_sequence_mutant_nc_m20.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m20.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m20_e1000_HGT_candidates.txt -t input_sequence_mutant_nc_m20.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m20_e2000_HGT_candidates.txt -t input_sequence_mutant_nc_m20.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m20_e5000_HGT_candidates.txt -t input_sequence_mutant_nc_m20.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m20_e10000_HGT_candidates.txt -t input_sequence_mutant_nc_m20.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m20_e50000_HGT_candidates.txt -t input_sequence_mutant_nc_m20.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m20_e1000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m20.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m20_e2000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m20.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m20_e5000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m20.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m20_e10000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m20.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m20_e50000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m20.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn


cd /srv/scratch/z5039045/MetaCHIP/m25
mkdir 7_get_gene_correlations_9_million
cp 6_MetaCHIP_9_million/output_ip90%_al200bp_c70%_e500bp/HGT_candidates_with_direction.txt 7_get_gene_correlations_9_million
cp 6_MetaCHIP_9_million/output_ip90_al200bp_c70_e1000bp/HGT_candidates_ET_validated.txt 7_get_gene_correlations_9_million
cp 6_MetaCHIP_9_million/combined.ffn 7_get_gene_correlations_9_million
cp ../assemblies/input_sequence_mutant_nc_m25.fasta 7_get_gene_correlations_9_million/
cd 7_get_gene_correlations_9_million
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs HGT_candidates_with_direction.txt -t input_sequence_mutant_nc_m25.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m25.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m25_e1000_HGT_candidates.txt -t input_sequence_mutant_nc_m25.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m25_e2000_HGT_candidates.txt -t input_sequence_mutant_nc_m25.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m25_e5000_HGT_candidates.txt -t input_sequence_mutant_nc_m25.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m25_e10000_HGT_candidates.txt -t input_sequence_mutant_nc_m25.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m25_e50000_HGT_candidates.txt -t input_sequence_mutant_nc_m25.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m25_e1000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m25.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m25_e2000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m25.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m25_e5000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m25.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m25_e10000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m25.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m25_e50000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m25.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn

	
cd /srv/scratch/z5039045/MetaCHIP/m30
mkdir 7_get_gene_correlations_9_million
cp 6_MetaCHIP_9_million/output_ip90%_al200bp_c70%_e500bp/HGT_candidates_with_direction.txt 7_get_gene_correlations_9_million
cp 6_MetaCHIP_9_million/output_ip90_al200bp_c70_e1000bp/HGT_candidates_ET_validated.txt 7_get_gene_correlations_9_million
cp 6_MetaCHIP_9_million/combined.ffn 7_get_gene_correlations_9_million
cp ../assemblies/input_sequence_mutant_nc_m30.fasta 7_get_gene_correlations_9_million/
cd 7_get_gene_correlations_9_million
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs HGT_candidates_with_direction.txt -t input_sequence_mutant_nc_m30.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m30.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m30_e1000_HGT_candidates.txt -t input_sequence_mutant_nc_m30.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m30_e2000_HGT_candidates.txt -t input_sequence_mutant_nc_m30.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m30_e5000_HGT_candidates.txt -t input_sequence_mutant_nc_m30.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m30_e10000_HGT_candidates.txt -t input_sequence_mutant_nc_m30.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m30_e50000_HGT_candidates.txt -t input_sequence_mutant_nc_m30.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m30_e1000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m30.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m30_e2000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m30.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m30_e5000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m30.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m30_e10000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m30.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m30_e50000_HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m30.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn

