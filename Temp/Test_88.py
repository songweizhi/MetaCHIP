
for each in open('/Users/songweizhi/Desktop/666.txt'):
    each = each.strip()
    #print('cd /srv/scratch/z5039045/MetaCHIP/NorthSea_samples/47bins_wd/prokka_output/%s' % each)
    #print('python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in %s.faa -t P &' % each)
    #print('')
    print('cp prokka_output/%s/%s_COG_results/func_stats.txt ./COG_annotation/%s_func_stats.txt &' % (each, each, each))



# python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in input.fasta -t P
