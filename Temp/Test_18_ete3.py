from ete3 import Tree

blast_results = '/Users/songweizhi/Desktop/assemblies/Ctg2ref_blast_results_m25.tab'
output =   open('/Users/songweizhi/Desktop/assemblies/Ctg2ref_m25.tab', 'w')

for each in open(blast_results):
    #print(each)
    each_split = each.strip().split()
    iden = float(each_split[2])
    align_len = int(each_split[3])
    query_len = int((each_split[12]))
    query = each_split[0]
    subject = each_split[1]
    #print(query_len)
    #print(iden)
    query_coverage = (align_len)/(query_len)

    if (iden >= 99) and (align_len >= 1000) and (query_coverage >= 0.5):
        #print(each.strip())
        #print(query_coverage)
        #print('%s\t%s' % (query, subject))
        output.write('%s\t%s\n' % (query, subject))
