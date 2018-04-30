import os
import shutil


# need to be changed
wd = '/Users/songweizhi/Desktop/wd'
bin_folder = '185bins_all'
checkm_wd = 'checkm_wd_universal'
bin_prefix = 'Refined_'
bin_ext = 'fasta'
completeness_cutoff = 40
contamination_cutoff = 0


# define file name
os.chdir(wd)
combined_qualities = 'combined_quality_%s.txt' % bin_folder
combined_qualities_filtered = 'combined_quality_%s_filtered.txt' % bin_folder
combined_qualities_qualified = 'combined_quality_%s_qualified_%s_%s.txt' % (bin_folder, completeness_cutoff, contamination_cutoff)
qualified_bin_folder = '%s_qualified_%s_%s' % (bin_folder, completeness_cutoff, contamination_cutoff)

# combine quality files
os.system('cat %s/%s/*.txt > %s.txt' % (bin_folder, checkm_wd, combined_qualities))


# create output folder
if os.path.isdir(qualified_bin_folder):
    shutil.rmtree(qualified_bin_folder)
    os.mkdir(qualified_bin_folder)
else:
    os.mkdir(qualified_bin_folder)


out_filtred = open(combined_qualities_filtered, 'w')
out_qualified = open(combined_qualities_qualified, 'w')
out_filtred.write('Bin_ID\tCompleteness\tContamination\n')
out_qualified.write('Bin_ID\tCompleteness\tContamination\n')
for each in open(combined_qualities):
    if each.startswith('  %s' % bin_prefix):
        each_split = each.strip().split(' ')
        each_split_new = []
        for each_element in each_split:
            if each_element != '':
                each_split_new.append(each_element)
        print(each_split_new)
        bin_id = each_split_new[0]
        completeness = each_split_new[12]
        contamination = each_split_new[13]
        out_filtered_content = '%s\t%s\t%s\n' % (bin_id, completeness, contamination)
        out_filtred.write(out_filtered_content)
        if (float(completeness) >= completeness_cutoff) and (float(contamination) <= contamination_cutoff):
            out_qualified.write(out_filtered_content)
            os.system('cp %s/%s.%s %s' % (bin_folder, bin_id, bin_ext, qualified_bin_folder))
out_filtred.close()
out_qualified.close()

