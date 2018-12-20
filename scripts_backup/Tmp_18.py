
for each in open('/Users/songweizhi/Desktop/combined_quality.txt'):


    if not (each.startswith('-')):
        if not (each.startswith('  Bin')):
            print(each.strip())

