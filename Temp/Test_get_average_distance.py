
n = 0
num_zero = 0
sum_all = 0
for each in open('/Users/songweizhi/Desktop/HG_NS/Human_gut_bins_distance_matrix/M.csv'):
    print(each)
    if not each.startswith(','):
        each_split = each.strip().split(',')
        bin_id = each_split[0]
        distance_list = each_split[1:]
        float_list = []
        for each_dis in distance_list:
            float_list.append(float(each_dis))
        print(float_list)
        n += len(float_list)
        sum_all += sum(float_list)
        num_zero += 1

print(n)
print(sum_all)
print(num_zero)
print(sum_all/(n - num_zero))
