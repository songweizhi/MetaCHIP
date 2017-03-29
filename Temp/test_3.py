

matches = open('/Users/songweizhi/Desktop/new_core/new_core.tab')

for each in matches:
    each_split = each.strip().split('\t')
    query = each_split[0]
    subject = each_split[1]
    if query != subject:
        print(each.strip())


