

for each in open('/Users/songweizhi/Desktop/test.txt'):
    each_split = each.strip().split('/')
    input_fasta = '/srv/scratch/z5039045/HgtSIM/reads_length/%s/combined.fasta' % (each.strip())
    #print(input_fasta)
    print('cd /srv/scratch/z5039045/HgtSIM/reads_length/%s' % (each.strip()))
    print('idba_ud --pre_correction --num_threads 1 --mink 20 --maxk 124 --step 20 --read %s --out %s_%s_k20-124' % (input_fasta, each_split[0], each_split[1]))



