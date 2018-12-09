import os


def unique_list_elements(list_input):

    list_output = []
    for each_element in list_input:
        if each_element not in list_output:
            list_output.append(each_element)

    return list_output


#release_version = 'r83'
release_version = 'r86'
wd = '/Users/songweizhi/Desktop/KelpBins/Taxon_Rank_Sankey/%s' % release_version
GTDB_output = '%s/gtdbtk.bac120.classification_528_%s.tsv' % (wd, release_version)

output_file = '%s_Class_to_Order.csv' % (release_version)
#output_file = '%s_Order_to_Family.csv' % (release_version)
#output_file = '%s_Class_to_Family.csv' % (release_version)


# path to get_sankey_plot.R
get_sankey_plot_Rscript = '~/PycharmProjects/Binning_refiner/get_sankey_plot.R'

# forward to wd
os.chdir(wd)

store_list = []
for each in open(GTDB_output):
    each_split = each.strip().split('\t')
    taxon_asign_list = each_split[1].split(';')

    assignment_full = []
    if len(each_split) == 1:
        assignment_full = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']

    elif (len(each_split) > 1) and (';' in each_split[1]):
        assignment = each_split[1].split(';')
        if len(assignment) == 7:
            assignment_full = assignment
        if len(assignment) == 6:
            assignment_full = assignment + ['s__']
        if len(assignment) == 5:
            assignment_full = assignment + ['g__', 's__']
        if len(assignment) == 4:
            assignment_full = assignment + ['f__', 'g__', 's__']
        if len(assignment) == 3:
            assignment_full = assignment + ['o__', 'f__', 'g__', 's__']
        if len(assignment) == 2:
            assignment_full = assignment + ['c__', 'o__', 'f__', 'g__', 's__']

    elif (len(each_split) > 1) and (';' not in each_split[1]):
        assignment_full = [each_split[1]] + ['p__', 'c__', 'o__', 'f__', 'g__', 's__']

    d = assignment_full[0]
    p = assignment_full[1]
    c = assignment_full[2]
    o = assignment_full[3]
    f = assignment_full[4]

    #store_list.append('%s,%s' % (p, c))
    store_list.append('%s,%s' % (c, o))
    #store_list.append('%s,%s' % (o, f))


store_list_uniq = unique_list_elements(store_list)
store_list_uniq_count_dict = {}
for each_key in store_list_uniq:
    store_list_uniq_count_dict[each_key] = store_list.count(each_key)


output_file_handle = open(output_file, 'w')
output_file_handle.write('C1,C2,Number\n')
for each in store_list_uniq_count_dict:
    for_out = '%s,%s\n' % (each, store_list_uniq_count_dict[each])
    #if not '__,' in for_out:
    output_file_handle.write(for_out)
output_file_handle.close()


# run Rscript
os.system('Rscript %s -f %s -y 1200' % (get_sankey_plot_Rscript, output_file))


# remove tmp files
os.system('rm %s' % output_file)
