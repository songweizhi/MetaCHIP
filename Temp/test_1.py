import os
import glob

os.chdir('/Users/songweizhi/Desktop')

files = '/Users/songweizhi/Desktop/all/*.txt'
file_list = [os.path.basename(file_name) for file_name in glob.glob(files)]
print(file_list)

print(len(file_list))


for each in file_list:
    each_fileted_handle = open('all_no_break/%s' % each, 'w')
    for each_line in open('all/%s' % each):
        print(each_line)
        each_line_split = each_line.strip().split('\t')
        end_break = each_line_split[5]
        print(end_break)
        if not end_break == 'yes':
            each_fileted_handle.write(each_line)

    each_fileted_handle.close()

