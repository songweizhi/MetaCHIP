
pwd_phylo_hmm_stat_txt = '/Users/songweizhi/Desktop/combined.txt'

for each in open(pwd_phylo_hmm_stat_txt):

    if not each.startswith('-'):
        if not each.startswith('  Bin'):
            print(each.strip())

