
for each in open('/Users/songweizhi/Desktop/combined_quality.txt'):
    if each.startswith('  NorthSea_bin'):
        print(each.strip())