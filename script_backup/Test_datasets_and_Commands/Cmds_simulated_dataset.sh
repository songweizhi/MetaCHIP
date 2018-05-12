
cd /Users/songweizhi/Desktop/MetaCHIP_wd
python3 /Users/songweizhi/PycharmProjects/MetaCHIP/MetaCHIP.py -cfg 


# test_1
cd /Users/songweizhi/Desktop/MetaCHIP_wd/168bins
python3 /Users/songweizhi/PycharmProjects/MetaCHIP/MetaCHIP.py -cfg config.txt



# optimize calculation step to 300 bp

# fq2fa
module load idba/1.1.3
fq2fa --merge combined_R1.fastq combined_R2.fastq combined.fasta

# get total length of simulated reads
module load python/3.5.2
python3 /srv/scratch/z5039045/Scripts/get_total_length.py -seq combined.fasta


########################### Performance on simulated datasets ############################

# simulate gene transfer with HgtSIM
python3 HgtSIM.py -t sequences_of_gene_transfers.fasta -i 0 -d distribution_of_transfers.txt -f selected_10_Betaproteobacteria -r 1-0-1-1 -x fna -lf TAGATGAGTGATTAGTTAGTTA -rf TAGATGAGTGATTAGTTAGTTA
python3 HgtSIM.py -t sequences_of_gene_transfers.fasta -i 5 -d distribution_of_transfers.txt -f selected_10_Betaproteobacteria -r 1-0-1-1 -x fna -lf TAGATGAGTGATTAGTTAGTTA -rf TAGATGAGTGATTAGTTAGTTA
python3 HgtSIM.py -t sequences_of_gene_transfers.fasta -i 10 -d distribution_of_transfers.txt -f selected_10_Betaproteobacteria -r 1-0-1-1 -x fna -lf TAGATGAGTGATTAGTTAGTTA -rf TAGATGAGTGATTAGTTAGTTA
python3 HgtSIM.py -t sequences_of_gene_transfers.fasta -i 15 -d distribution_of_transfers.txt -f selected_10_Betaproteobacteria -r 1-0-1-1 -x fna -lf TAGATGAGTGATTAGTTAGTTA -rf TAGATGAGTGATTAGTTAGTTA
python3 HgtSIM.py -t sequences_of_gene_transfers.fasta -i 20 -d distribution_of_transfers.txt -f selected_10_Betaproteobacteria -r 1-0-1-1 -x fna -lf TAGATGAGTGATTAGTTAGTTA -rf TAGATGAGTGATTAGTTAGTTA
python3 HgtSIM.py -t sequences_of_gene_transfers.fasta -i 25 -d distribution_of_transfers.txt -f selected_10_Betaproteobacteria -r 1-0-1-1 -x fna -lf TAGATGAGTGATTAGTTAGTTA -rf TAGATGAGTGATTAGTTAGTTA
python3 HgtSIM.py -t sequences_of_gene_transfers.fasta -i 30 -d distribution_of_transfers.txt -f selected_10_Betaproteobacteria -r 1-0-1-1 -x fna -lf TAGATGAGTGATTAGTTAGTTA -rf TAGATGAGTGATTAGTTAGTTA

###### positive control ######
# Run Prokka 
module load perl/5.20.1
module load infernal/1.1.1
module load blast+/2.2.31
module load hmmer/3.1b2
module load prodigal/2.6.3
module load tbl2asn/25.3
module load parallel/20160222
module load prokka/1.12
prokka --force --prefix BAD --locustag BAD --strain BAD --outdir BAD BAD.fna
# run all_vs_all blast
module load blast+/2.6.0
makeblastdb -in combined.ffn -dbtype nucl -parse_seqids
blastn -query combined.ffn -db combined.ffn -out all_vs_all_ffn.tab -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn
# run MetaCHIP 
module load python/3.5.2
module load blast+/2.6.0
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/MetaCHIP.py -cfg config.txt


################################ optimize sequencing depth ###############################

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

# get assembler recovered HGTs
module load python/3.5.2
module load blast+/2.6.0

# m0
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers_iden100.py -a m0_IDBA_UD_3_million_k20-124.fasta -t input_sequence_mutant_nc_m0.fasta -minf 1000 -d distribution_of_transfers.txt -combined_genomes combined_reference_genomes_m0.fna
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers_iden100.py -a m0_IDBA_UD_6_million_k20-124.fasta -t input_sequence_mutant_nc_m0.fasta -minf 1000 -d distribution_of_transfers.txt -combined_genomes combined_reference_genomes_m0.fna
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers_iden100.py -a m0_IDBA_UD_9_million_k20-124.fasta -t input_sequence_mutant_nc_m0.fasta -minf 1000 -d distribution_of_transfers.txt -combined_genomes combined_reference_genomes_m0.fna
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers_iden100.py -a m0_IDBA_UD_12_million_k20-124.fasta -t input_sequence_mutant_nc_m0.fasta -minf 1000 -d distribution_of_transfers.txt -combined_genomes combined_reference_genomes_m0.fna
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers_iden100.py -a m0_metaSPAdes_3_million_scaffolds.fasta -t input_sequence_mutant_nc_m0.fasta -minf 1000 -d distribution_of_transfers.txt -combined_genomes combined_reference_genomes_m0.fna
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers_iden100.py -a m0_metaSPAdes_6_million_scaffolds.fasta -t input_sequence_mutant_nc_m0.fasta -minf 1000 -d distribution_of_transfers.txt -combined_genomes combined_reference_genomes_m0.fna
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers_iden100.py -a m0_metaSPAdes_9_million_scaffolds.fasta -t input_sequence_mutant_nc_m0.fasta -minf 1000 -d distribution_of_transfers.txt -combined_genomes combined_reference_genomes_m0.fna
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers_iden100.py -a m0_metaSPAdes_12_million_scaffolds.fasta -t input_sequence_mutant_nc_m0.fasta -minf 1000 -d distribution_of_transfers.txt -combined_genomes combined_reference_genomes_m0.fna
# m5
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m5_IDBA_UD_3_million_k20-124.fasta -t input_sequence_mutant_nc_m5.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m5_IDBA_UD_6_million_k20-124.fasta -t input_sequence_mutant_nc_m5.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m5_IDBA_UD_9_million_k20-124.fasta -t input_sequence_mutant_nc_m5.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m5_IDBA_UD_12_million_k20-124.fasta -t input_sequence_mutant_nc_m5.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m5_metaSPAdes_3_million_scaffolds.fasta -t input_sequence_mutant_nc_m5.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m5_metaSPAdes_6_million_scaffolds.fasta -t input_sequence_mutant_nc_m5.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m5_metaSPAdes_9_million_scaffolds.fasta -t input_sequence_mutant_nc_m5.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m5_metaSPAdes_12_million_scaffolds.fasta -t input_sequence_mutant_nc_m5.fasta -minf 1000
# m10
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m10_IDBA_UD_3_million_k20-124.fasta -t input_sequence_mutant_nc_m10.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m10_IDBA_UD_6_million_k20-124.fasta -t input_sequence_mutant_nc_m10.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m10_IDBA_UD_9_million_k20-124.fasta -t input_sequence_mutant_nc_m10.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m10_IDBA_UD_12_million_k20-124.fasta -t input_sequence_mutant_nc_m10.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m10_metaSPAdes_3_million_scaffolds.fasta -t input_sequence_mutant_nc_m10.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m10_metaSPAdes_6_million_scaffolds.fasta -t input_sequence_mutant_nc_m10.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m10_metaSPAdes_9_million_scaffolds.fasta -t input_sequence_mutant_nc_m10.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m10_metaSPAdes_12_million_scaffolds.fasta -t input_sequence_mutant_nc_m10.fasta -minf 1000
# m15
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m15_IDBA_UD_3_million_k20-124.fasta -t input_sequence_mutant_nc_m15.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m15_IDBA_UD_6_million_k20-124.fasta -t input_sequence_mutant_nc_m15.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m15_IDBA_UD_9_milliona_k20-124.fasta -t input_sequence_mutant_nc_m15.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m15_IDBA_UD_12_million_k20-124.fasta -t input_sequence_mutant_nc_m15.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m15_metaSPAdes_3_million_scaffolds.fasta -t input_sequence_mutant_nc_m15.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m15_metaSPAdes_6_million_scaffolds.fasta -t input_sequence_mutant_nc_m15.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m15_metaSPAdes_9_million_scaffolds.fasta -t input_sequence_mutant_nc_m15.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m15_metaSPAdes_12_million_scaffolds.fasta -t input_sequence_mutant_nc_m15.fasta -minf 1000
# m20
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m20_IDBA_UD_3_million_k20-124.fasta -t input_sequence_mutant_nc_m20.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m20_IDBA_UD_6_million_k20-124.fasta -t input_sequence_mutant_nc_m20.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m20_IDBA_UD_9_million_k20-124.fasta -t input_sequence_mutant_nc_m20.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m20_IDBA_UD_12_million_k20-124.fasta -t input_sequence_mutant_nc_m20.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m20_metaSPAdes_3_million_scaffolds.fasta -t input_sequence_mutant_nc_m20.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m20_metaSPAdes_6_million_scaffolds.fasta -t input_sequence_mutant_nc_m20.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m20_metaSPAdes_9_million_scaffolds.fasta -t input_sequence_mutant_nc_m20.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m20_metaSPAdes_12_million_scaffolds.fasta -t input_sequence_mutant_nc_m20.fasta -minf 1000
# m25
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m25_IDBA_UD_3_million_k20-124.fasta -t input_sequence_mutant_nc_m25.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m25_IDBA_UD_6_million_k20-124.fasta -t input_sequence_mutant_nc_m25.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m25_IDBA_UD_9_million_k20-124.fasta -t input_sequence_mutant_nc_m25.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m25_IDBA_UD_12_million_k20-124.fasta -t input_sequence_mutant_nc_m25.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m25_metaSPAdes_3_million_scaffolds.fasta -t input_sequence_mutant_nc_m25.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m25_metaSPAdes_6_million_scaffolds.fasta -t input_sequence_mutant_nc_m25.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m25_metaSPAdes_9_million_scaffolds.fasta -t input_sequence_mutant_nc_m25.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m25_metaSPAdes_12_million_scaffolds.fasta -t input_sequence_mutant_nc_m25.fasta -minf 1000
# m30
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m30_IDBA_UD_3_million_k20-124.fasta -t input_sequence_mutant_nc_m30.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m30_IDBA_UD_6_million_k20-124.fasta -t input_sequence_mutant_nc_m30.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m30_IDBA_UD_9_million_k20-124.fasta -t input_sequence_mutant_nc_m30.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m30_IDBA_UD_12_million_k20-124.fasta -t input_sequence_mutant_nc_m30.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m30_metaSPAdes_3_million_scaffolds.fasta -t input_sequence_mutant_nc_m30.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m30_metaSPAdes_6_million_scaffolds.fasta -t input_sequence_mutant_nc_m30.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m30_metaSPAdes_9_million_scaffolds.fasta -t input_sequence_mutant_nc_m30.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a m30_metaSPAdes_12_million_scaffolds.fasta -t input_sequence_mutant_nc_m30.fasta -minf 1000
# mNC
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a mNC_IDBA_UD_9_million_k20-124.fasta -t input_sequence_mutant_nc_mNC.fasta -minf 1000
python3 /srv/scratch/z5039045/HgtSIM/Assembler_recovered_transfers.py -a mNC_metaSPAdes_9_million_scaffolds.fasta -t input_sequence_mutant_nc_mNC.fasta -minf 1000


# summary
Group	Reads_num		metaSPAdes	IDBA_UD 
m0		3 million		51			24
m0		6 million		71			11
m0		9 million		57			3
m0		12 million		51			2

m5		3 million		3			11
m5		6 million		5			19
m5		9 million		8			23
m5		12 million		11			30

m10		3 million		6			17
m10		6 million		22			91
m10		9 million		34			98
m10		12 million		35			98

m15		3 million		10			32
m15		6 million		63			96
m15		9 million		90			98
m15		12 million		89			99

m20		3 million		19			42
m20		6 million		88			98
m20		9 million		98			100
m20		12 million		99			100

m25		3 million		29			39
m25		6 million		92			99
m25		9 million		99			100
m25		12 million		99			100

m30		3 million		28			41
m30		6 million		98			99
m30		9 million		100			100
m30		12 million		98			100


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
cd /srv/scratch/z5039045/MetaCHIP/m0
cd /srv/scratch/z5039045/MetaCHIP/m5
cd /srv/scratch/z5039045/MetaCHIP/m10

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
# module load blast+/2.6.0
# module load mafft/7.310
# module load hmmer/3.1b2
# module load fasttree/2.1.7
module load python/2.7.12
module load ete3/3.0.0b36
virtualenv --system-site-packages mypythonenv
. mypythonenv/bin/activate
export PATH=~/anaconda_ete/bin:$PATH;

# m0
cd /srv/scratch/z5039045/MetaCHIP/m0/6_MetaCHIP_9_million
python /srv/scratch/z5039045/Softwares/MetaCHIP/MetaCHIP.py -cfg config.txt
python /srv/scratch/z5039045/Softwares/MetaCHIP/Explicit_Tree.py -cfg config.txt
# m5
cd /srv/scratch/z5039045/MetaCHIP/m5/6_MetaCHIP_9_million
python /srv/scratch/z5039045/Softwares/MetaCHIP/MetaCHIP.py -cfg config.txt
python /srv/scratch/z5039045/Softwares/MetaCHIP/Explicit_Tree.py -cfg config.txt
# m10
cd /srv/scratch/z5039045/MetaCHIP/m10/6_MetaCHIP_9_million
python /srv/scratch/z5039045/Softwares/MetaCHIP/MetaCHIP.py -cfg config.txt
python /srv/scratch/z5039045/Softwares/MetaCHIP/Explicit_Tree.py -cfg config.txt
# m15
cd /srv/scratch/z5039045/MetaCHIP/m15/6_MetaCHIP_9_million
python /srv/scratch/z5039045/Softwares/MetaCHIP/MetaCHIP.py -cfg config.txt
python /srv/scratch/z5039045/Softwares/MetaCHIP/Explicit_Tree.py -cfg config.txt
# m20
cd /srv/scratch/z5039045/MetaCHIP/m20/6_MetaCHIP_9_million
python /srv/scratch/z5039045/Softwares/MetaCHIP/MetaCHIP.py -cfg config.txt
python /srv/scratch/z5039045/Softwares/MetaCHIP/Explicit_Tree.py -cfg config.txt
# m25
cd /srv/scratch/z5039045/MetaCHIP/m25/6_MetaCHIP_9_million
python /srv/scratch/z5039045/Softwares/MetaCHIP/MetaCHIP.py -cfg config.txt
python /srv/scratch/z5039045/Softwares/MetaCHIP/Explicit_Tree.py -cfg config.txt
# m30
cd /srv/scratch/z5039045/MetaCHIP/m30/6_MetaCHIP_9_million
python /srv/scratch/z5039045/Softwares/MetaCHIP/MetaCHIP.py -cfg config.txt
python /srv/scratch/z5039045/Softwares/MetaCHIP/Explicit_Tree.py -cfg config.txt


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





cd /srv/scratch/z5039045/MetaCHIP/m0/positive_control_m0
python /srv/scratch/z5039045/Softwares/MetaCHIP/Explicit_Tree.py -cfg config.txt
cd /srv/scratch/z5039045/MetaCHIP/m5/positive_control_m5
python /srv/scratch/z5039045/Softwares/MetaCHIP/Explicit_Tree.py -cfg config.txt
cd /srv/scratch/z5039045/MetaCHIP/m10/positive_control_m10
python /srv/scratch/z5039045/Softwares/MetaCHIP/Explicit_Tree.py -cfg config.txt
cd /srv/scratch/z5039045/MetaCHIP/m15/positive_control_m15
python /srv/scratch/z5039045/Softwares/MetaCHIP/Explicit_Tree.py -cfg config.txt
cd /srv/scratch/z5039045/MetaCHIP/m20/positive_control_m20
python /srv/scratch/z5039045/Softwares/MetaCHIP/Explicit_Tree.py -cfg config.txt
cd /srv/scratch/z5039045/MetaCHIP/m25/positive_control_m25
python /srv/scratch/z5039045/Softwares/MetaCHIP/Explicit_Tree.py -cfg config.txt
cd /srv/scratch/z5039045/MetaCHIP/m30/positive_control_m30
python /srv/scratch/z5039045/Softwares/MetaCHIP/Explicit_Tree.py -cfg config.txt









# summary
Group	Reads_num		Size			Depth		metaSPAdes	IDBA_UD MetaCHIP(SPAdes)	MetaCHIP(idba)
m0		3 million		429.41 Mbp		5.67x		51			24		8
m0		4.2 million		601.3 Mbp		7.93x		57			26		34
m0		4.5 million		644.28 Mbp		8.50x		66			24		31
m0		6 million		859.09 Mbp		11.33x		71			11		40
m0		7.5 million		1073.65 Mbp		14.16x		60			4		37
m0		9 million		1288.53 Mbp		17.0x		57			3		27
m0		10.5 million	1503.09	Mbp		19.83x		67			1		42
m0		12 million		1717.76	Mbp		22.66x		51			2		39

m5		3 million									3			11
m5		4.2 million									9			10
m5		4.5 million									12			13
m5		6 million									5			19
m5		7.5 million									12			23
m5		9 million									8			23
m5		10.5 million								15			25
m5		12 million									11			30
m5		18 million									11			38
m5		24 million									3			34
m5		30 million									0			36

m10		3 million									6			17
m10		4.2 million									11			49
m10		4.5 million									12			60
m10		6 million									22			91
m10		7.5 million									36			93
m10		9 million									34			98
m10		10.5 million								42			97
m10		12 million									35			98

m15		3 million									10			32
m15		6 million									63			96
m15		9 million									90			98
m15		12 million									89			99
m15		18 million									89			98

m20		3 million									19			42
m20		6 million									88			98
m20		9 million									98			100
m20		12 million									99			100
m20		18 million									98			99

m25		3 million									29			39
m25		6 million									92			99
m25		9 million									99			100
m25		12 million									99			100
m25		18 million									99			100

m30		3 million									28			41
m30		6 million									98			99
m30		9 million									100			100
m30		12 million									98			100
m30		18 million									98			99



# get COG annotation of 20 genomes
module unload python/2.7.12
module load python/3.5.2
module load perl/5.20.1
module load blast+/2.6.0
cd /srv/scratch/z5039045/MetaCHIP/mNC/positive_control_mNC/prokka_output/AAM
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in AAM.faa -t P &
cd /srv/scratch/z5039045/MetaCHIP/mNC/positive_control_mNC/prokka_output/AKV
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in AKV.faa -t P &
cd /srv/scratch/z5039045/MetaCHIP/mNC/positive_control_mNC/prokka_output/AMAC
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in AMAC.faa -t P &
cd /srv/scratch/z5039045/MetaCHIP/mNC/positive_control_mNC/prokka_output/AMAU
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in AMAU.faa -t P &
cd /srv/scratch/z5039045/MetaCHIP/mNC/positive_control_mNC/prokka_output/AMS
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in AMS.faa -t P &
cd /srv/scratch/z5039045/MetaCHIP/mNC/positive_control_mNC/prokka_output/ARL
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in ARL.faa -t P &
cd /srv/scratch/z5039045/MetaCHIP/mNC/positive_control_mNC/prokka_output/ARS
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in ARS.faa -t P &
cd /srv/scratch/z5039045/MetaCHIP/mNC/positive_control_mNC/prokka_output/ASJ
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in ASJ.faa -t P &
cd /srv/scratch/z5039045/MetaCHIP/mNC/positive_control_mNC/prokka_output/ASN
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in ASN.faa -t P &
cd /srv/scratch/z5039045/MetaCHIP/mNC/positive_control_mNC/prokka_output/ATM
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in ATM.faa -t P &
cd /srv/scratch/z5039045/MetaCHIP/mNC/positive_control_mNC/prokka_output/BAD
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in BAD.faa -t P &
cd /srv/scratch/z5039045/MetaCHIP/mNC/positive_control_mNC/prokka_output/BDS
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in BDS.faa -t P &
cd /srv/scratch/z5039045/MetaCHIP/mNC/positive_control_mNC/prokka_output/BGC
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in BGC.faa -t P &
cd /srv/scratch/z5039045/MetaCHIP/mNC/positive_control_mNC/prokka_output/BHS
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in BHS.faa -t P &
cd /srv/scratch/z5039045/MetaCHIP/mNC/positive_control_mNC/prokka_output/BNM
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in BNM.faa -t P &
cd /srv/scratch/z5039045/MetaCHIP/mNC/positive_control_mNC/prokka_output/BRT
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in BRT.faa -t P &
cd /srv/scratch/z5039045/MetaCHIP/mNC/positive_control_mNC/prokka_output/BSA
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in BSA.faa -t P &
cd /srv/scratch/z5039045/MetaCHIP/mNC/positive_control_mNC/prokka_output/BSD
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in BSD.faa -t P &
cd /srv/scratch/z5039045/MetaCHIP/mNC/positive_control_mNC/prokka_output/BSL
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in BSL.faa -t P &
cd /srv/scratch/z5039045/MetaCHIP/mNC/positive_control_mNC/prokka_output/BTK
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in BTK.faa -t P &



cd /srv/scratch/z5039045/MetaCHIP/mNC/positive_control_mNC
mkdir COG_annotation
cp prokka_output/AAM/AAM_COG_results/func_stats.txt ./COG_annotation/AAM_func_stats.txt
cp prokka_output/AKV/AKV_COG_results/func_stats.txt ./COG_annotation/AKV_func_stats.txt
cp prokka_output/AMAC/AMAC_COG_results/func_stats.txt ./COG_annotation/AMAC_func_stats.txt
cp prokka_output/AMAU/AMAU_COG_results/func_stats.txt ./COG_annotation/AMAU_func_stats.txt
cp prokka_output/AMS/AMS_COG_results/func_stats.txt ./COG_annotation/AMS_func_stats.txt
cp prokka_output/ARL/ARL_COG_results/func_stats.txt ./COG_annotation/ARL_func_stats.txt
cp prokka_output/ARS/ARS_COG_results/func_stats.txt ./COG_annotation/ARS_func_stats.txt
cp prokka_output/ASJ/ASJ_COG_results/func_stats.txt ./COG_annotation/ASJ_func_stats.txt
cp prokka_output/ASN/ASN_COG_results/func_stats.txt ./COG_annotation/ASN_func_stats.txt
cp prokka_output/ATM/ATM_COG_results/func_stats.txt ./COG_annotation/ATM_func_stats.txt
cp prokka_output/BAD/BAD_COG_results/func_stats.txt ./COG_annotation/BAD_func_stats.txt
cp prokka_output/BDS/BDS_COG_results/func_stats.txt ./COG_annotation/BDS_func_stats.txt
cp prokka_output/BGC/BGC_COG_results/func_stats.txt ./COG_annotation/BGC_func_stats.txt
cp prokka_output/BHS/BHS_COG_results/func_stats.txt ./COG_annotation/BHS_func_stats.txt
cp prokka_output/BNM/BNM_COG_results/func_stats.txt ./COG_annotation/BNM_func_stats.txt
cp prokka_output/BRT/BRT_COG_results/func_stats.txt ./COG_annotation/BRT_func_stats.txt
cp prokka_output/BSA/BSA_COG_results/func_stats.txt ./COG_annotation/BSA_func_stats.txt
cp prokka_output/BSD/BSD_COG_results/func_stats.txt ./COG_annotation/BSD_func_stats.txt
cp prokka_output/BSL/BSL_COG_results/func_stats.txt ./COG_annotation/BSL_func_stats.txt
cp prokka_output/BTK/BTK_COG_results/func_stats.txt ./COG_annotation/BTK_func_stats.txt


module load python/3.5.2
module load perl/5.20.1
module load blast+/2.6.0
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in HGT_candidates.fasta -t N
python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in HGT_candidates_ET_validated.fasta -t N







# get the number of recovered gene transfers by assemblers from dornor, recipient or both
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




python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/Assembler_recovered_transfers_iden100.py -a m0_metaSPAdes_9_million_scaffolds.fasta -t input_sequence_mutant_nc_m0.fasta -minf 1000 -d distribution_of_transfers.txt -combined_genomes combined_reference_genomes_m0.fna 



python3 /srv/scratch/z5039045/Binning_refiner/CheckM_runner.py -1 Renamed_refined_bins_m20 -2 Renamed_refined_bins_m25 -r Renamed_refined_bins_m30 -pbs -qsub





########## assess SCG tree ##########

# run Prokka
module unload python/2.7.12
module load python/3.5.2
python3 /srv/scratch/z5039045/Softwares/MetaCHIP/InPrep1_Prokka.py -genome_folder genomes -genome_extension fna

# get species tree and grouping input genomes
module unload python/3.5.2 
module load python/2.7.12
module load blast+/2.6.0
module load mafft/7.310
module load hmmer/3.1b2
module load fasttree/2.1.7
python /srv/scratch/z5039045/Softwares/MetaCHIP/Get_species_tree.py -prokka_output prokka_output -hmm /srv/scratch/z5039045/Softwares/MetaCHIP/phylo.hmm 
python /srv/scratch/z5039045/Softwares/MetaCHIP/Get_species_tree.py -prokka_output prokka_output -hmm /srv/scratch/z5039045/Softwares/MetaCHIP/PCG74.hmm 
python /srv/scratch/z5039045/Softwares/MetaCHIP/Get_species_tree.py -prokka_output prokka_output -hmm /srv/scratch/z5039045/Softwares/MetaCHIP/PCG82.hmm 



########## determine end break length ##########

# copy files into working directory
cp m0/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e1000bp/HGT_candidates.fasta determine_end_break_length/m0_e1000_HGT_candidates.fasta &
cp m0/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e2000bp/HGT_candidates.fasta determine_end_break_length/m0_e2000_HGT_candidates.fasta &
cp m0/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e5000bp/HGT_candidates.fasta determine_end_break_length/m0_e5000_HGT_candidates.fasta &
cp m0/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e10000bp/HGT_candidates.fasta determine_end_break_length/m0_e10000_HGT_candidates.fasta &
cp m0/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e50000bp/HGT_candidates.fasta determine_end_break_length/m0_e50000_HGT_candidates.fasta &
cp m5/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e1000bp/HGT_candidates.fasta determine_end_break_length/m5_e1000_HGT_candidates.fasta &
cp m5/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e2000bp/HGT_candidates.fasta determine_end_break_length/m5_e2000_HGT_candidates.fasta &
cp m5/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e5000bp/HGT_candidates.fasta determine_end_break_length/m5_e5000_HGT_candidates.fasta &
cp m5/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e10000bp/HGT_candidates.fasta determine_end_break_length/m5_e10000_HGT_candidates.fasta &
cp m5/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e50000bp/HGT_candidates.fasta determine_end_break_length/m5_e50000_HGT_candidates.fasta &
cp m10/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e1000bp/HGT_candidates.fasta determine_end_break_length/m10_e1000_HGT_candidates.fasta &
cp m10/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e2000bp/HGT_candidates.fasta determine_end_break_length/m10_e2000_HGT_candidates.fasta &
cp m10/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e5000bp/HGT_candidates.fasta determine_end_break_length/m10_e5000_HGT_candidates.fasta &
cp m10/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e10000bp/HGT_candidates.fasta determine_end_break_length/m10_e10000_HGT_candidates.fasta &
cp m10/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e50000bp/HGT_candidates.fasta determine_end_break_length/m10_e50000_HGT_candidates.fasta &
cp m15/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e1000bp/HGT_candidates.fasta determine_end_break_length/m15_e1000_HGT_candidates.fasta &
cp m15/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e2000bp/HGT_candidates.fasta determine_end_break_length/m15_e2000_HGT_candidates.fasta &
cp m15/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e5000bp/HGT_candidates.fasta determine_end_break_length/m15_e5000_HGT_candidates.fasta &
cp m15/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e10000bp/HGT_candidates.fasta determine_end_break_length/m15_e10000_HGT_candidates.fasta &
cp m15/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e50000bp/HGT_candidates.fasta determine_end_break_length/m15_e50000_HGT_candidates.fasta &
cp m20/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e1000bp/HGT_candidates.fasta determine_end_break_length/m20_e1000_HGT_candidates.fasta &
cp m20/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e2000bp/HGT_candidates.fasta determine_end_break_length/m20_e2000_HGT_candidates.fasta &
cp m20/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e5000bp/HGT_candidates.fasta determine_end_break_length/m20_e5000_HGT_candidates.fasta &
cp m20/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e10000bp/HGT_candidates.fasta determine_end_break_length/m20_e10000_HGT_candidates.fasta &
cp m20/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e50000bp/HGT_candidates.fasta determine_end_break_length/m20_e50000_HGT_candidates.fasta &
cp m25/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e1000bp/HGT_candidates.fasta determine_end_break_length/m25_e1000_HGT_candidates.fasta &
cp m25/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e2000bp/HGT_candidates.fasta determine_end_break_length/m25_e2000_HGT_candidates.fasta &
cp m25/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e5000bp/HGT_candidates.fasta determine_end_break_length/m25_e5000_HGT_candidates.fasta &
cp m25/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e10000bp/HGT_candidates.fasta determine_end_break_length/m25_e10000_HGT_candidates.fasta &
cp m25/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e50000bp/HGT_candidates.fasta determine_end_break_length/m25_e50000_HGT_candidates.fasta &
cp m30/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e1000bp/HGT_candidates.fasta determine_end_break_length/m30_e1000_HGT_candidates.fasta &
cp m30/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e2000bp/HGT_candidates.fasta determine_end_break_length/m30_e2000_HGT_candidates.fasta &
cp m30/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e5000bp/HGT_candidates.fasta determine_end_break_length/m30_e5000_HGT_candidates.fasta &
cp m30/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e10000bp/HGT_candidates.fasta determine_end_break_length/m30_e10000_HGT_candidates.fasta &
cp m30/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e50000bp/HGT_candidates.fasta determine_end_break_length/m30_e50000_HGT_candidates.fasta &
cp m0/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e1000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m0_e1000_HGT_candidates_ET_validated.fasta &
cp m0/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e2000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m0_e2000_HGT_candidates_ET_validated.fasta &
cp m0/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e5000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m0_e5000_HGT_candidates_ET_validated.fasta &
cp m0/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e10000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m0_e10000_HGT_candidates_ET_validated.fasta &
cp m0/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e50000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m0_e50000_HGT_candidates_ET_validated.fasta &
cp m5/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e1000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m5_e1000_HGT_candidates_ET_validated.fasta &
cp m5/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e2000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m5_e2000_HGT_candidates_ET_validated.fasta &
cp m5/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e5000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m5_e5000_HGT_candidates_ET_validated.fasta &
cp m5/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e10000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m5_e10000_HGT_candidates_ET_validated.fasta &
cp m5/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e50000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m5_e50000_HGT_candidates_ET_validated.fasta &
cp m10/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e1000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m10_e1000_HGT_candidates_ET_validated.fasta &
cp m10/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e2000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m10_e2000_HGT_candidates_ET_validated.fasta &
cp m10/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e5000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m10_e5000_HGT_candidates_ET_validated.fasta &
cp m10/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e10000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m10_e10000_HGT_candidates_ET_validated.fasta &
cp m10/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e50000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m10_e50000_HGT_candidates_ET_validated.fasta &
cp m15/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e1000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m15_e1000_HGT_candidates_ET_validated.fasta &
cp m15/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e2000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m15_e2000_HGT_candidates_ET_validated.fasta &
cp m15/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e5000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m15_e5000_HGT_candidates_ET_validated.fasta &
cp m15/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e10000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m15_e10000_HGT_candidates_ET_validated.fasta &
cp m15/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e50000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m15_e50000_HGT_candidates_ET_validated.fasta &
cp m20/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e1000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m20_e1000_HGT_candidates_ET_validated.fasta &
cp m20/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e2000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m20_e2000_HGT_candidates_ET_validated.fasta &
cp m20/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e5000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m20_e5000_HGT_candidates_ET_validated.fasta &
cp m20/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e10000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m20_e10000_HGT_candidates_ET_validated.fasta &
cp m20/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e50000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m20_e50000_HGT_candidates_ET_validated.fasta &
cp m25/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e1000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m25_e1000_HGT_candidates_ET_validated.fasta &
cp m25/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e2000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m25_e2000_HGT_candidates_ET_validated.fasta &
cp m25/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e5000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m25_e5000_HGT_candidates_ET_validated.fasta &
cp m25/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e10000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m25_e10000_HGT_candidates_ET_validated.fasta &
cp m25/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e50000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m25_e50000_HGT_candidates_ET_validated.fasta &
cp m30/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e1000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m30_e1000_HGT_candidates_ET_validated.fasta &
cp m30/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e2000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m30_e2000_HGT_candidates_ET_validated.fasta &
cp m30/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e5000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m30_e5000_HGT_candidates_ET_validated.fasta &
cp m30/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e10000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m30_e10000_HGT_candidates_ET_validated.fasta &
cp m30/6_MetaCHIP_9_million/output_ip90_al200bp_c70_e50000bp/HGT_candidates_ET_validated.fasta determine_end_break_length/m30_e50000_HGT_candidates_ET_validated.fasta &

# get results
python3 Step_6_2_get_gene_correlations_iden100.py -predicted_HGTs m0_e1000_HGT_candidates.fasta -t input_sequence_mutant_nc_m0.fasta -d distribution_of_transfers.txt
python3 Step_6_2_get_gene_correlations_iden100.py -predicted_HGTs m0_e1000_HGT_candidates_ET_validated.fasta -t input_sequence_mutant_nc_m0.fasta -d distribution_of_transfers.txt
python3 Step_6_2_get_gene_correlations_iden100.py -predicted_HGTs m0_e2000_HGT_candidates.fasta -t input_sequence_mutant_nc_m0.fasta -d distribution_of_transfers.txt
python3 Step_6_2_get_gene_correlations_iden100.py -predicted_HGTs m0_e2000_HGT_candidates_ET_validated.fasta -t input_sequence_mutant_nc_m0.fasta -d distribution_of_transfers.txt
python3 Step_6_2_get_gene_correlations_iden100.py -predicted_HGTs m0_e5000_HGT_candidates.fasta -t input_sequence_mutant_nc_m0.fasta -d distribution_of_transfers.txt
python3 Step_6_2_get_gene_correlations_iden100.py -predicted_HGTs m0_e5000_HGT_candidates_ET_validated.fasta -t input_sequence_mutant_nc_m0.fasta -d distribution_of_transfers.txt
python3 Step_6_2_get_gene_correlations_iden100.py -predicted_HGTs m0_e10000_HGT_candidates.fasta -t input_sequence_mutant_nc_m0.fasta -d distribution_of_transfers.txt
python3 Step_6_2_get_gene_correlations_iden100.py -predicted_HGTs m0_e10000_HGT_candidates_ET_validated.fasta -t input_sequence_mutant_nc_m0.fasta -d distribution_of_transfers.txt
python3 Step_6_2_get_gene_correlations_iden100.py -predicted_HGTs m0_e50000_HGT_candidates.fasta -t input_sequence_mutant_nc_m0.fasta -d distribution_of_transfers.txt
python3 Step_6_2_get_gene_correlations_iden100.py -predicted_HGTs m0_e50000_HGT_candidates_ET_validated.fasta -t input_sequence_mutant_nc_m0.fasta -d distribution_of_transfers.txt




python3 Step_6_2_get_gene_correlations.py -predicted_HGTs m10_e1000_HGT_candidates.fasta -t input_sequence_mutant_nc_m10.fasta -d distribution_of_transfers.txt
python3 Step_6_2_get_gene_correlations.py -predicted_HGTs m10_e1000_HGT_candidates_ET_validated.fasta -t input_sequence_mutant_nc_m10.fasta -d distribution_of_transfers.txt






python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs m10_e1000_HGT_candidates.fasta -t input_sequence_mutant_nc_m20.fasta -d distribution_of_transfers.txt
python3 /srv/scratch/z5039045/MetaCHIP/Step_6_2_get_gene_correlations.py -predicted_HGTs HGT_candidates_ET_validated.txt -t input_sequence_mutant_nc_m20.fasta -d distribution_of_transfers.txt -combined_ffn combined.ffn
python3 Step_6_2_get_gene_correlations.py -predicted_HGTs m10_e10000_HGT_candidates.fasta -t input_sequence_mutant_nc_m10.fasta -d distribution_of_transfers.txt
python3 Step_6_2_get_gene_correlations.py -predicted_HGTs m10_e10000_HGT_candidates_ET_validated.fasta -t input_sequence_mutant_nc_m10.fasta -d distribution_of_transfers.txt





