import os

mutation_levels = [0, 5, 10, 15, 20, 25, 30]

# wd = '/Users/songweizhi/Desktop/HgtSIM_wd'
# os.chdir(wd)

bootstrap_num = 1
while bootstrap_num <= 10:
    for mutation_level in mutation_levels:
        cmd_1_1 = 'mkdir MetaCHIP_wd/bootstrap%s/bootstrap%s_m%s' % (bootstrap_num, bootstrap_num, mutation_level)
        cmd_1_2 = 'mkdir MetaCHIP_wd/bootstrap%s/bootstrap%s_m%s/imput_genome_m%s' % (bootstrap_num, bootstrap_num, mutation_level, mutation_level)
        cmd_2 = 'cp HgtSIM_wd/bootstrap%s/bootstrap%s_outputs_%s_1-0-1-1/Genomes_with_transfers/*.fna MetaCHIP_wd/bootstrap%s/bootstrap%s_m%s/imput_genome_m%s/' % (bootstrap_num, bootstrap_num, mutation_level, bootstrap_num, bootstrap_num, mutation_level, mutation_level)
        cmd_3_between_genus = 'cp genome_folder_SB/*.fna MetaCHIP_wd/bootstrap%s/bootstrap%s_m%s/imput_genome_m%s/' % (bootstrap_num, bootstrap_num, mutation_level, mutation_level)
        #cmd_3_between_class = 'cp genome_folder_alpha/*.fna MetaCHIP_wd/bootstrap%s/bootstrap%s_m%s/imput_genome_m%s/' % (bootstrap_num, bootstrap_num, mutation_level, mutation_level)
        cmd_4 = 'cp GTDB_results.txt MetaCHIP_wd/bootstrap%s/bootstrap%s_m%s/' % (bootstrap_num, bootstrap_num, mutation_level)
        print(cmd_1_1)
        print(cmd_1_2)
        print(cmd_2)
        print(cmd_3_between_genus)
        #print(cmd_3_between_class)
        print(cmd_4)
        #os.system(cmd_1_1)
        #os.system(cmd_1_2)
        #os.system(cmd_2)
        #os.system(cmd_3_between_genus)
        #os.system(cmd_3_between_class)
        #os.system(cmd_4)
        print()

    bootstrap_num += 1
