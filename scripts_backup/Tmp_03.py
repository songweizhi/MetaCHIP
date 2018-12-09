import matplotlib.pyplot as plt
import numpy as np

summary_file = '/Users/songweizhi/Desktop/summary_class_level.txt'
summary_file = '/Users/songweizhi/Desktop/summary_genus_level.txt'


recovered_BM_list = [[], [], [], [], [], [], []]
recovered_PG_list = [[], [], [], [], [], [], []]
recovered_RD_list = [[], [], [], [], [], [], []]
for each in open(summary_file):
    each_split = each.strip().split('\t')
    bootstrap = each_split[0].split('_')[0]
    mutation_level = each_split[0].split('_')[1]
    recovered_BM = int(each_split[1])
    recovered_PG = int(each_split[2])
    recovered_RD = int(each_split[3])

    if mutation_level == 'm0':
        recovered_BM_list[0].append(recovered_BM)
        recovered_PG_list[0].append(recovered_PG)
        recovered_RD_list[0].append(recovered_RD)

    if mutation_level == 'm5':
        recovered_BM_list[1].append(recovered_BM)
        recovered_PG_list[1].append(recovered_PG)
        recovered_RD_list[1].append(recovered_RD)

    if mutation_level == 'm10':
        recovered_BM_list[2].append(recovered_BM)
        recovered_PG_list[2].append(recovered_PG)
        recovered_RD_list[2].append(recovered_RD)

    if mutation_level == 'm15':
        recovered_BM_list[3].append(recovered_BM)
        recovered_PG_list[3].append(recovered_PG)
        recovered_RD_list[3].append(recovered_RD)

    if mutation_level == 'm20':
        recovered_BM_list[4].append(recovered_BM)
        recovered_PG_list[4].append(recovered_PG)
        recovered_RD_list[4].append(recovered_RD)

    if mutation_level == 'm25':
        recovered_BM_list[5].append(recovered_BM)
        recovered_PG_list[5].append(recovered_PG)
        recovered_RD_list[5].append(recovered_RD)

    if mutation_level == 'm30':
        recovered_BM_list[6].append(recovered_BM)
        recovered_PG_list[6].append(recovered_PG)
        recovered_RD_list[6].append(recovered_RD)


recovered_BM_matrix = np.asmatrix(recovered_BM_list)
recovered_PG_matrix = np.asmatrix(recovered_PG_list)
recovered_RD_matrix = np.asmatrix(recovered_RD_list)

recovered_BM_matrix_t = recovered_BM_matrix.transpose()
recovered_PG_matrix_t = recovered_PG_matrix.transpose()
recovered_RD_matrix_t = recovered_RD_matrix.transpose()

BM_col_mean = np.mean(recovered_BM_matrix_t, axis=0)
PG_col_mean = np.mean(recovered_PG_matrix_t, axis=0)
RD_col_mean = np.mean(recovered_RD_matrix_t, axis=0)

BM_col_std = np.std(recovered_BM_matrix_t, axis=0)
PG_col_std = np.std(recovered_PG_matrix_t, axis=0)
RD_col_std = np.std(recovered_RD_matrix_t, axis=0)

BM_col_mean_list = BM_col_mean.tolist()
PG_col_mean_list = PG_col_mean.tolist()
RD_col_mean_list = RD_col_mean.tolist()

BM_col_std_list = BM_col_std.tolist()
PG_col_std_list = PG_col_std.tolist()
RD_col_std_list = RD_col_std.tolist()


print(BM_col_mean_list[0])
print(PG_col_mean_list[0])
print(RD_col_mean_list[0])
print()

print(BM_col_std_list[0])
print(PG_col_std_list[0])
print(RD_col_std_list[0])
print()




# plt.plot(BM_col_mean_list[0])
# #plt.errorbar(7, BM_col_std_list[0])
#
# plt.plot(PG_col_mean_list[0])
# plt.plot(RD_col_mean_list[0])
#
#
# plt.ylabel('some numbers')
# plt.show()
#
