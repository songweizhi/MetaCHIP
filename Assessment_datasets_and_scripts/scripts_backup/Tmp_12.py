import os
import glob

wd = '/Users/songweizhi/Desktop/555555/'

BM_output_folder = 'BM_HGTs'
BM_output_folder = 'MetaBAT_bins'
os.chdir(wd)


BM_output_file_re = '%s/*.txt' % BM_output_folder
BM_output_file_list = [os.path.basename(file_name) for file_name in glob.glob(BM_output_file_re)]


HGT_occurence_dict = {}
for BM_output_file in BM_output_file_list:
    pwd_BM_output_file = '%s/%s' % (BM_output_folder, BM_output_file)
    taxon_rank = BM_output_file.split('_')[1][0]
    for BM_HGT in open(pwd_BM_output_file):
        if not BM_HGT.startswith('Gene_1'):
            BM_HGT_split = BM_HGT.strip().split('\t')
            concatenated = '%s___%s' % (BM_HGT_split[0], BM_HGT_split[1])

            if concatenated not in HGT_occurence_dict:
                HGT_occurence_dict[concatenated] = [taxon_rank]
            else:
                HGT_occurence_dict[concatenated].append(taxon_rank)

n = 0
for each_hgt in HGT_occurence_dict:
    for_out = '%s\t\t%s' % (each_hgt, ''.join(HGT_occurence_dict[each_hgt]))
    print(for_out)
    n += len(HGT_occurence_dict[each_hgt])

print(n)


print(len(HGT_occurence_dict))


# NorthSea_bin029_00160___NorthSea_bin052_00351	o
# NorthSea_bin029_00160___NorthSea_bin067_02110	f

# NorthSea_bin029_01349___NorthSea_bin052_00753	o
# NorthSea_bin029_01349___NorthSea_bin090_00799	f

# NorthSea_bin029_02254___NorthSea_bin052_01493	o
# NorthSea_bin029_02254___NorthSea_bin067_01843	f

# NorthSea_bin040_00403___NorthSea_bin053_00702	c
# NorthSea_bin040_00403___NorthSea_bin083_01550	o

# NorthSea_bin050_00353___NorthSea_bin029_00160	g
# NorthSea_bin050_00353___NorthSea_bin052_00351	o
# NorthSea_bin050_00353___NorthSea_bin067_02110	f

# NorthSea_bin050_01455___NorthSea_bin029_00538	g
# NorthSea_bin050_01455___NorthSea_bin067_01091	f

# NorthSea_bin107_00955___NorthSea_bin028_00753	c
# NorthSea_bin107_00955___NorthSea_bin069_01663	o



