import os
import glob


wd = '/Users/songweizhi/Desktop/555555/TT_90MGs'
PG_output_folder = 'TT_90MGs_PG'
combined_output = 'TT_90MGs_PG.txt'


os.chdir(wd)


PG_output_file_re = '%s/*.txt' % PG_output_folder
PG_output_file_list_with_path = [('%s/%s' % (PG_output_folder, os.path.basename(file_name))) for file_name in glob.glob(PG_output_file_re)]


for each_file in PG_output_file_list_with_path:

    end_match_num = 0
    full_length_match_num = 0
    PG_validated_num = 0
    qualified_num = 0
    for each_hgt in open(each_file):
        if not each_hgt.startswith('Gene_1'):
            each_hgt_split = each_hgt.strip().split('\t')
            end_match = each_hgt_split[5]
            full_length_match = each_hgt_split[6]
            direction = each_hgt_split[7]
            if direction != 'NA':
                PG_validated_num += 1
                if end_match == 'yes':
                    end_match_num += 1
                if full_length_match == 'yes':
                    full_length_match_num += 1
                if (end_match == 'no') and (full_length_match == 'no'):
                    qualified_num += 1

    for_print = '%s\t%s\t%s\t%s\t%s' % (each_file, PG_validated_num, end_match_num, full_length_match_num, qualified_num)
    print(for_print)




end_match_num = 0
full_length_match_num = 0
PG_validated_num = 0
qualified_num = 0
for hgt in open(combined_output):
    if not hgt.startswith('Gene_1'):
        each_hgt_split = hgt.strip().split('\t')

        #print(each_hgt_split)

        end_match = each_hgt_split[4]
        full_length_match = each_hgt_split[5]
        direction = each_hgt_split[6]
        if direction != 'NA':
            PG_validated_num += 1
            if end_match == 'yes':
                end_match_num += 1
            if full_length_match == 'yes':
                full_length_match_num += 1
            if (end_match == 'no') and (full_length_match == 'no'):
                qualified_num += 1

for_print = '%s\t%s\t%s\t%s\t%s' % ('combined', PG_validated_num, end_match_num, full_length_match_num, qualified_num)
print(for_print)



