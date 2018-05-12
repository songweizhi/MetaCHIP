#  http://stackoverflow.com/questions/14391959/heatmap-in-matplotlib-with-pcolor

import string
import matplotlib.pyplot as plt
import numpy as np

group_idens = open("/Users/weizhisong/Desktop/working_directory/0_output_3rd_average.txt")
#out = open("/Users/weizhisong/Desktop/working_directory/for_plot_in_r.txt",'w')
group_numbers = 12


def get_group_identity_dic(input_file):
    group_identity_dic = {}
    for each_line in input_file:
        each_line_split = each_line.strip().split('\t')
        group_name = each_line_split[0]
        identity = each_line_split[1]
        group_identity_dic[group_name] = identity
    return group_identity_dic


def get_identity_list(group_iden_dic, group_number):
    all_alphabet_list = list(string.ascii_uppercase)
    sub_alphabet_list = all_alphabet_list[0: group_number]
    members = []
    members_iden_str = []
    members_iden_flo = []
    #out.write (',' + '%s' % ','.join(sub_alphabet_list) + '\n') # The first title line
    #out.write ('%s' % ','.join(sub_alphabet_list) + '\n') # The first title line
    for alphabet_r in sub_alphabet_list:
        r_list = []
        r_list_iden_str = []
        r_list_iden_flo = []
        for alphabet_c in sub_alphabet_list:
            combined = alphabet_r + '_' + alphabet_c
            r_list.append(combined)
            if combined in group_iden_dic:
                r_list_iden_str.append(group_iden_dic[combined])
                r_list_iden_flo.append(float(group_iden_dic[combined]))
            else:
                r_list_iden_str.append('0')
                r_list_iden_flo.append(0)
        members.append(r_list)
        members_iden_str.append(r_list_iden_str)
        members_iden_flo.append(r_list_iden_flo)
        #out.write (alphabet_r + ',' + '%s' % ','.join(r_list_iden_str) + '\n')
        #out.write ('%s' % ','.join(r_list_iden_str) + '\n')
    return members_iden_flo


def get_heatmap_plot(iden_list, group_number):
    data = np.array(iden_list)
    column_labels = list(string.ascii_uppercase)[0: group_number]
    row_labels = list(string.ascii_uppercase)[0: group_number]
    fig, ax = plt.subplots()
    heatmap = ax.pcolor(data)
    plt.colorbar(heatmap)
    ax.set_xticks(np.arange(data.shape[0])+0.5, minor=False) # put the major ticks at the middle of each cell
    ax.set_yticks(np.arange(data.shape[1])+0.5, minor=False) # put the major ticks at the middle of each cell
    ax.invert_yaxis() #setup the direction of Y-axis
    ax.xaxis.tick_top() #put X-axis on the top of the figure
    ax.set_xticklabels(row_labels, minor=False)
    ax.set_yticklabels(column_labels, minor=False)
    plt.show()


group_identity_dic = get_group_identity_dic(group_idens)
list_a = get_identity_list(group_identity_dic,group_numbers)
get_heatmap_plot(list_a,group_numbers)
