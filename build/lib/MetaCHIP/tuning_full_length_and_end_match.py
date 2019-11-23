

def check_match_direction(blast_hit_splitted):
    query_start = int(blast_hit_splitted[6])
    query_end = int(blast_hit_splitted[7])
    subject_start = int(blast_hit_splitted[8])
    subject_end = int(blast_hit_splitted[9])
    query_direction = query_end - query_start
    subject_direction = subject_end - subject_start

    same_match_direction = True
    if ((query_direction > 0) and (subject_direction < 0)) or ((query_direction < 0) and (subject_direction > 0)):
        same_match_direction = False

    return same_match_direction


def check_full_length_and_end_match(qualified_ctg_match_list, identity_cutoff):

    ######################################## check full length match ########################################

    query_len = int(qualified_ctg_match_list[0][12])
    subject_len = int(qualified_ctg_match_list[0][13])

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
            if query_matched[1] > current_query_end:
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
            if subject_matched[1] > current_subject_end:
                subject_matched_len_total += subject_matched[1] - current_subject_end

        elif subject_matched[0] == current_subject_end:
            subject_matched_len_total += subject_matched[1] - current_subject_end


    # get total coverage for query and subject
    query_cov_total = query_matched_len_total/query_len
    subject_cov_total = subject_matched_len_total/subject_len

    # get match category
    match_category = 'normal'
    best_hit_end_gap_len = 200
    gap_cutoff_for_concatenating = 300

    # full length match: coverage cutoff 90%
    if (query_cov_total >= 0.90) or (subject_cov_total >= 0.90):
        match_category = 'full_length_match'


    ######################################## check end match ########################################

    else:

        # read in best hit information
        best_hit = qualified_ctg_match_list[0]
        best_hit_identity = float(best_hit[2])
        best_hit_query_start = int(best_hit[6])
        best_hit_query_end = int(best_hit[7])
        best_hit_subject_start = int(best_hit[8])
        best_hit_subject_end = int(best_hit[9])
        query_len = int(best_hit[12])
        subject_len = int(best_hit[13])
        best_hit_same_direction = check_match_direction(best_hit)

        # concatenate continuously matched blocks with gap less than 200bp
        matched_block_query_start = best_hit_query_start
        matched_block_query_end = best_hit_query_end
        matched_block_subject_start = best_hit_subject_start
        matched_block_subject_end = best_hit_subject_end

        if best_hit_identity >= identity_cutoff:

            for matched_block in qualified_ctg_match_list[1:]:

                current_block_identity = float(matched_block[2])
                current_block_direction = check_match_direction(matched_block)

                # if identity difference <= 1 and has same match direction with the best hit
                if (-6 <= (best_hit_identity - current_block_identity) <= 6) and (current_block_direction == best_hit_same_direction):

                    current_query_start = int(matched_block[6])
                    current_query_end = int(matched_block[7])
                    current_subject_start = int(matched_block[8])
                    current_subject_end = int(matched_block[9])

                    if best_hit_same_direction is True:

                        # situation 1
                        if ((current_query_start >= matched_block_query_start) and (current_query_end <= matched_block_query_end)) and ((current_subject_start >= matched_block_subject_start) and (current_subject_end <= matched_block_subject_end)):
                            pass  # do nothing

                        # situation 2
                        if ((current_query_start > matched_block_query_start) and (current_query_end > matched_block_query_end)) and ((current_subject_start > matched_block_subject_start) and (current_subject_end > matched_block_subject_end)) and (-gap_cutoff_for_concatenating <= (current_query_start - matched_block_query_end) <= gap_cutoff_for_concatenating) and (-gap_cutoff_for_concatenating <= (current_subject_start - matched_block_subject_end) <= gap_cutoff_for_concatenating):
                            matched_block_query_end = current_query_end
                            matched_block_subject_end = current_subject_end

                        # situation 3
                        if ((current_query_start < matched_block_query_start) and (current_query_end < matched_block_query_end)) and ((current_subject_start < matched_block_subject_start) and (current_subject_end < matched_block_subject_end)) and (-gap_cutoff_for_concatenating <= (matched_block_query_start - current_query_end) <= gap_cutoff_for_concatenating) and (-gap_cutoff_for_concatenating <= (matched_block_subject_start - current_subject_end) <= gap_cutoff_for_concatenating):
                            matched_block_query_start = current_query_start
                            matched_block_subject_start = current_subject_start

                    if best_hit_same_direction is False:

                        # situation 1
                        if ((current_query_start >= matched_block_query_start) and (current_query_end <= matched_block_query_end)) and ((current_subject_start <= matched_block_subject_start) and (current_subject_end >= matched_block_subject_end)):
                            pass

                        # situation 2
                        if ((current_query_start > matched_block_query_start) and (current_query_end > matched_block_query_end)) and ((current_subject_start < matched_block_subject_start) and (current_subject_end < matched_block_subject_end)) and (-gap_cutoff_for_concatenating <= (current_query_start - matched_block_query_end) <= gap_cutoff_for_concatenating) and (-gap_cutoff_for_concatenating <= (matched_block_subject_end - current_subject_start) <= gap_cutoff_for_concatenating):
                            matched_block_query_end = current_query_end
                            matched_block_subject_end = current_subject_end

                        # situation 3
                        if ((current_query_start < matched_block_query_start) and (current_query_end < matched_block_query_end)) and ((current_subject_start > matched_block_subject_start) and (current_subject_end > matched_block_subject_end)) and (-gap_cutoff_for_concatenating <= (matched_block_query_start - current_query_end) <= gap_cutoff_for_concatenating) and (-gap_cutoff_for_concatenating <= (current_subject_end - matched_block_subject_start) <= gap_cutoff_for_concatenating):
                            matched_block_query_start = current_query_start
                            matched_block_subject_start = current_subject_start


            ######################################## check end_match ########################################

            # situation 1
            if (best_hit_same_direction is True) and (query_len - matched_block_query_end <= best_hit_end_gap_len) and (matched_block_subject_start <= best_hit_end_gap_len):
                match_category = 'end_match'

            # situation 2
            elif (best_hit_same_direction is True) and (matched_block_query_start <= best_hit_end_gap_len) and (subject_len - matched_block_subject_end <= best_hit_end_gap_len):
                match_category = 'end_match'

            # situation 3
            elif (best_hit_same_direction is False) and (query_len - matched_block_query_end <= best_hit_end_gap_len) and (subject_len - matched_block_subject_start <= best_hit_end_gap_len):
                match_category = 'end_match'

            # situation 4
            elif (best_hit_same_direction is False) and (matched_block_query_start <= best_hit_end_gap_len) and (matched_block_subject_end <= best_hit_end_gap_len):
                match_category = 'end_match'

    return match_category


#output_c_full_len = '/Users/songweizhi/Desktop/000/RHO2_bin_73_00772___RHO3_bin_59_01813/RHO2_bin_73_00772___RHO3_bin_59_01813_full_length.txt'

#output_c_full_len = '/Users/songweizhi/Desktop/000/COS2_bin_21_00105___APA_bin_43_02125/COS2_bin_21_00105___APA_bin_43_02125_full_length.txt'

output_c_full_len = '/Users/songweizhi/Desktop/000/COS3_bin_21_01438___IRC_PAM_SB0661_bin_58_00227/COS3_bin_21_01438___IRC_PAM_SB0661_bin_58_00227_full_length.txt'

#output_c_full_len = '/Users/songweizhi/Desktop/000/COS4_bin_26_02662___IRC_PAM_SB0668_bin_42_00545/COS4_bin_26_02662___IRC_PAM_SB0668_bin_42_00545_full_length.txt'



# get qualified_ctg_match_list
min_ctg_match_aln_len = 100
qualified_ctg_match_list = []
for blast_hit in open(output_c_full_len):
    blast_hit_split = blast_hit.strip().split('\t')
    align_len = int(blast_hit_split[3])
    if align_len >= min_ctg_match_aln_len:
        qualified_ctg_match_list.append(blast_hit_split)

if len(qualified_ctg_match_list) == 0:
    match_category = 'normal'
else:
    match_category = check_full_length_and_end_match(qualified_ctg_match_list, 90)


print(match_category)
