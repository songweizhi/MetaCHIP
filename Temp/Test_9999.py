import os

# get id to category dict
id2cate_dict = {}
id2fun_dict = {}

for each in open('/Users/songweizhi/Desktop/ARDBflatFiles/cog.tab'):
    cog_id = 'COG%s' % each.strip().split('\t')[0]
    cog_category = each.strip().split('\t')[1]
    cog_fun = each.strip().split('\t')[2]
    id2cate_dict[cog_id] = cog_category
    id2fun_dict[cog_id] = cog_fun


ar_cog_list = []
for each2 in open('/Users/songweizhi/Desktop/ARDBflatFiles/type2cog.tab'):
    cog_id2 = 'COG%s' % each2.strip().split('\t')[1]
    if cog_id2 not in ar_cog_list:
        ar_cog_list.append(cog_id2)

out = open('/Users/songweizhi/Desktop/AR_COG_category.txt', 'w')
for each3 in ar_cog_list:
    out.write('%s\t%s\t%s\n' % (id2cate_dict[each3], each3, id2fun_dict[each3]))
out.close()

os.system('cat /Users/songweizhi/Desktop/AR_COG_category.txt | sort > /Users/songweizhi/Desktop/AR_COG_category_sorted.txt')


for each4 in ['COG0254', 'COG0256', 'COG0291', 'COG0292', 'COG0261', 'COG0222', 'COG0228']:
    print(id2cate_dict[each4])

