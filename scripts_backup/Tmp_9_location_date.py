
genome_id_file = '/Users/songweizhi/Desktop/KelpBins/genome_id_list.txt'


# get 528 genome is list
genome_528_list = []
for genome in open(genome_id_file):
    genome_528_list.append(genome.strip())

date_list = ['050416', '060416', '080416',  # 2016-04: 05, 06, 08
             '070616', '140616',            # 2016-06: 07, 14
             '050716', '070716',            # 2016-07: 05, 07
             '110816', '170816',            # 2016-08: 11, 17
             '061016',                      # 2016-10: 06
             '141216', '161216',            # 2016-12: 14, 16
             '070217', '080217',            # 2017-02: 07, 08
             '040417', '050417',            # 2017-04: 04, 05
             '050617', '130617']            # 2017-06: 05, 13

location_list = ['BH', 'BI', 'CB', 'SH']
host = 'ER'

# get matrix
print('Date\tBH\tBI\tCB\tSH')
for date in date_list:

    current_date_stats = [date]
    for location in location_list:
        tag = '%s_%s_%s' % (location, host, date)
        current_tag_genome_num = 0
        for genome in genome_528_list:
            if tag in genome:
                current_tag_genome_num += 1
        current_date_stats.append(current_tag_genome_num)

    current_date_stats_str = [str(i) for i in current_date_stats]
    for_write = '\t'.join(current_date_stats_str)
    print(for_write)
