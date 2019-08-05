

file_in = '/Users/songweizhi/Desktop/MetaBAT138bins_o29_HGTs_BM.txt'
#file_in = '/Users/songweizhi/Desktop/MetaBAT138bins_o29_HGTs_PG_validated.txt'


file_in = '/Users/songweizhi/Desktop/NorthSea37bins_o16_HGTs_BM.txt'
#file_in = '/Users/songweizhi/Desktop/NorthSea37bins_o16_HGTs_PG_validated.txt'




num_0 = 0
num_5 = 0
num_10 = 0
num_15 = 0
num_20 = 0
num_25 = 0
num_30 = 0
for each in open(file_in):
    if not each.startswith('Gene_1'):


        identity = float(each.strip().split('\t')[4])

        gv = 100 - identity

        print(gv)

        if gv <= 2.5:
            num_0 += 1
        elif 2.5 < gv <= 7.5:
            num_5 += 1
        elif 7.5 < gv <= 12.5:
            num_10 += 1
        elif 12.5 < gv <= 17.5:
            num_15 += 1
        elif 17.5 < gv <= 22.5:
            num_20 += 1
        elif 22.5 < gv <= 27.5:
            num_25 += 1
        elif 27.5 < gv <= 32.5:
            num_30 += 1




print('%s %s %s %s %s %s %s' % (num_0, num_5, num_10, num_15, num_20, num_25, num_30))




