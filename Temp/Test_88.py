
for each in open('/Users/songweizhi/Desktop/26bins.txt'):
    each = each.strip()
    print('cp ./prokka_output/%s/%s_COG_results/func_stats.txt ./COG_annotation_bins/%s_func_stats.txt' % (each, each, each))






# python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in input.fasta -t P
