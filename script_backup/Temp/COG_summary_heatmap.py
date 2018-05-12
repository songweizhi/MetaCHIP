import os
import glob
import pandas
import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# python3 ~/PycharmProjects/MetaCHIP/Temp/COG_summary_heatmap.py -in func_stats_files -out PBTR_ML.jpg -in_percent

def turn_to_percentage(number_list):
    number_list_percent = []
    for each_element in number_list:
        each_element_percent = float("{0:.2f}".format(each_element / sum(number_list)))
        number_list_percent.append(each_element_percent)
    return number_list_percent


parser = argparse.ArgumentParser()

parser.add_argument('-in',
                    required=True,
                    help='folder holds the func_stats.txt file for each genome')

parser.add_argument('-out',
                    required=True,
                    help='output image')

parser.add_argument('-columns',
                    required=False,
                    help='the order of columns')

parser.add_argument('-in_percent',
                    required=False,
                    action="store_true",
                    help='in percent')

args = vars(parser.parse_args())
input_folder = args['in']
image_title = args['out']
column_order = args['columns']
in_percent = args['in_percent']


arrange_order_column = []
if column_order != None:
    arrange_order_column = column_order.split(',')

arrange_order_row = ['J', 'A', 'K', 'L', 'B', 'D', 'Y', 'V', 'T', 'M', 'N', 'Z', 'W', 'U', 'O', 'C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q', 'R', 'S']
arrange_order_row_reversed = arrange_order_row[::-1]


input_files = '%s/*.txt' % input_folder
file_list = [os.path.basename(file_name) for file_name in glob.glob(input_files)]
file_list = sorted(file_list)


category_id_list = sorted(arrange_order_row)
genome_name_list = []
category_num_list = []
for each_file in file_list:
    genome_name = '_'.join(each_file.split('_')[:-2])
    genome_name_list.append(genome_name)
    current_category_num_list = []
    for each_category in open('%s/%s' % (input_folder, each_file)):
        each_category_count = int(each_category.strip().split('\t')[2])
        current_category_num_list.append(each_category_count)

    # turn absolute number to percentage if specified
    if in_percent == 1:
        current_category_num_list = turn_to_percentage(current_category_num_list)

    category_num_list.append(current_category_num_list)
category_num_arrary = np.transpose(np.array(category_num_list))


# check the number of input genomes and provided column order
if column_order != None:
    if len(genome_name_list) != len(arrange_order_column):
        print('Dected %s input genomes and specified %s columns, The input genomes and specified columns are inconsistent.' % (len(genome_name_list), len(arrange_order_column)))
        exit()


# add row and column name to dataframe
category_num_df = pandas.DataFrame(category_num_arrary, index=category_id_list, columns=genome_name_list)


# re-index row order
category_num_df = category_num_df.reindex(arrange_order_row_reversed)

# re-index column order
genome_name_list_sorted = sorted(genome_name_list)
if column_order == None:
    category_num_df = category_num_df.reindex(genome_name_list_sorted, axis=1)
else:
    category_num_df = category_num_df.reindex(arrange_order_column, axis=1)

print(category_num_df)

# get the minimal value in the dataframe
min_value = 0
n = 0
for each_min in category_num_df.min():
    if n == 0:
        min_value = each_min
    else:
        if each_min < min_value:
            min_value = each_min
    n += 1

# get the maximal value in the dataframe
max_value = 0
m = 0
for each_max in category_num_df.max():
    if m == 0:
        max_value = each_max
    else:
        if each_max > max_value:
            max_value = each_max
    m += 1


# set color bar
cdict = {'red':  ((0.0, 0.0, 0.0),   # no red at 0
                  (1.0, 1.0, 1.0)),  # set red to 1
        'green': ((0.0, 0.0, 0.0),   # no green at 0
                  (1.0, 0.0, 0.0)),  # no green at 1
        'blue':  ((0.0, 0.0, 0.0),   # no blue at 0
                  (1.0, 0.0, 0.0))   # no blue at 1
       }
GnRd = colors.LinearSegmentedColormap('GnRd', cdict)
fig, ax = plt.subplots(1)
p = ax.pcolormesh(category_num_df, cmap=GnRd, vmin=min_value, vmax=max_value)
fig.colorbar(p,ax=ax, shrink=1) # orientation="horizontal"

# set x and y ticks
ax.set_yticks(np.arange(category_num_df.shape[0])+0.5, minor=False)  # direction = 'out'

ax.set_xticks(np.arange(category_num_df.shape[1])+0.5, minor=False)

# set x and y labels
ax.set_yticklabels(arrange_order_row_reversed, minor=False)

# set the direction of ticks
ax.tick_params(direction='out', length=0)

if column_order == None:
    ax.set_xticklabels(genome_name_list_sorted, minor=False, rotation=90)
else:
    ax.set_xticklabels(arrange_order_column, minor=False, rotation=90)

ax.xaxis.tick_top()
ax.yaxis.tick_right()

# set figure size
row_num = category_num_df.shape[0]
column_num = category_num_df.shape[1]
#fig.set_size_inches(column_num + 3, row_num)

plt.savefig('%s.jpg' % (image_title), dpi = 600)
plt.close()
