

min_ctg_match_aln_len = 100

blast_results = '/Users/songweizhi/Desktop/check_end/rh_hh_Refined_15_00163___rh_lh_Refined_25_02576/rh_hh_Refined_15_00163___rh_lh_Refined_25_02576_full_length.txt'

qualified_ctg_match_list = []
for blast_hit in open(blast_results):
    blast_hit_split = blast_hit.strip().split('\t')
    align_len = int(blast_hit_split[3])
    query_len = int(blast_hit_split[12])
    subject_len = int(blast_hit_split[13])
    if align_len >= min_ctg_match_aln_len:
        qualified_ctg_match_list.append(blast_hit_split)


def check_full_lenght_and_end_match(qualified_ctg_match_list, identity_cutoff):

    # get position list of matched regions
    query_matched_region_list = []
    subject_matched_region_list = []
    for ctg_match in qualified_ctg_match_list:
        query_start = int(ctg_match[6])
        query_end = int(ctg_match[7])
        subject_start = int(ctg_match[8])
        subject_end = int(ctg_match[9])
        query_matched_region = sorted([query_start, query_end])
        subject_matched_region = sorted([subject_start, subject_end])
        query_matched_region_list.append(query_matched_region)
        subject_matched_region_list.append(subject_matched_region)

    # get total length of query matched regions
    query_matched_len_total = 0
    current_query_end = 0
    for query_matched in sorted(query_matched_region_list):
        if query_matched_len_total == 0:
            query_matched_len_total = query_matched[1] - query_matched[0] + 1
            current_query_end = query_matched[1]
        elif query_matched[0] > current_query_end:
            query_matched_len_total += query_matched[1] - query_matched[0] + 1
        elif query_matched[0] < current_query_end:
            query_matched_len_total += query_matched[1] - current_query_end
        elif query_matched[0] == current_query_end:
            query_matched_len_total += query_matched[1] - current_query_end

    # get total length of subject matched regions
    subject_matched_len_total = 0
    current_subject_end = 0
    for subject_matched in sorted(subject_matched_region_list):
        if subject_matched_len_total == 0:
            subject_matched_len_total = subject_matched[1] - subject_matched[0] + 1
            current_subject_end = subject_matched[1]
        elif subject_matched[0] > current_subject_end:
            subject_matched_len_total += subject_matched[1] - subject_matched[0] + 1
        elif subject_matched[0] < current_subject_end:
            subject_matched_len_total += subject_matched[1] - current_subject_end
        elif subject_matched[0] == current_subject_end:
            subject_matched_len_total += subject_matched[1] - current_subject_end

    # get total coverage for query and subject
    query_cov_total = query_matched_len_total/query_len
    subject_cov_total = subject_matched_len_total/subject_len

    # get match category
    match_category = 'normal'
    best_hit_end_gap_len = 200

    # full length match: coverage cutoff 95%
    if (query_cov_total >= 0.95) or (subject_cov_total >= 0.95):
        match_category = 'full_length_match'

    else:
        best_hit = qualified_ctg_match_list[0]
        best_hit_identity = float(best_hit[2])

        if best_hit_identity >= identity_cutoff:
            best_hit_query_start =   int(best_hit[6])
            best_hit_query_end =     int(best_hit[7])
            best_hit_subject_start = int(best_hit[8])
            best_hit_subject_end =   int(best_hit[9])
            best_hit_query_len =     int(best_hit[12])
            best_hit_subject_len =   int(best_hit[13])

            # check best_hit match direction
            best_hit_query_direction = best_hit_query_end - best_hit_query_start
            best_hit_subject_direction = best_hit_subject_end - best_hit_subject_start
            best_hit_same_match_direction = True
            if ((best_hit_query_direction > 0) and (best_hit_subject_direction < 0)) or ((best_hit_query_direction < 0) and (best_hit_subject_direction > 0)):
                best_hit_same_match_direction = False

            # situation 1
            if (best_hit_same_match_direction is True) and (best_hit_query_len - best_hit_query_end <= best_hit_end_gap_len) and (best_hit_subject_start <= best_hit_end_gap_len):
                match_category = 'end_match'

            # situation 2
            elif (best_hit_same_match_direction is True) and (best_hit_query_start <= best_hit_end_gap_len) and (best_hit_subject_len - best_hit_subject_end <= best_hit_end_gap_len):
                match_category = 'end_match'

            # situation 3
            elif (best_hit_same_match_direction is False) and (best_hit_query_len - best_hit_query_end <= best_hit_end_gap_len) and (best_hit_subject_len - best_hit_subject_start <= best_hit_end_gap_len):
                match_category = 'end_match'

            # situation 4
            elif (best_hit_same_match_direction is False) and (best_hit_query_start <= best_hit_end_gap_len) and (best_hit_subject_end <= best_hit_end_gap_len):
                match_category = 'end_match'

    return match_category

qualified_ctg_match_list = []
identity_cutoff = 99
match_category = check_full_lenght_and_end_match(qualified_ctg_match_list, identity_cutoff)
print(match_category)

