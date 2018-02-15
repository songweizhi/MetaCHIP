import os

os.chdir('/Users/songweizhi/Desktop/test_end_break')


for each_match in open('all_vs_all_500bp.tab'):
    match_split = each_match.strip().split('\t')
    query = match_split[0]
    subject = match_split[1]
    identity = float(match_split[2])
    align_len = int(match_split[3])
    query_len = int(match_split[12])
    subject_len = int(match_split[13])
    if (identity >= 98) and (align_len >= 400) and (query != subject):
        print(each_match)
















