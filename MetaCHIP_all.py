import sys
import argparse
from scripts.PI import PI
from scripts.BM import BM
from scripts_backup.PG_all import PG

to_do = '''

add option -force
remove outline from Soil199_grouping_o46.png 
flk plot, remove matches shorter than 50bp
not print if disabled: [2018-12-09 21:54:32] Plotting flanking regions with 16 cores
ignore GTDB results header line
report progress for every 100 processes (if n/100 is int)
move config dict to MetaCHIP.py
move uclust step to PG.py
error ocurred during plot if no HGT validated by PG
run mafft in fast mode if sequence number higher than 200 
change text direction for circos plot
report getting genome size

'''


def print_main_help():

    help_message = ''' 
             ...::: MetaCHIP :::...
        
    HGT detection modules:
       PI               -> Prepare Input files 
       BM               -> Best-Match approach 
       PG               -> PhyloGenetic approach
    
    Plot modules:
       plot_tree        -> Plot newick tree
       plot_taxon       -> show taxon correlations with sankey plot
      
    Other modules:
       parallel_blastn  -> run all-vs-all blastn in parallel
       get_gene_cluster -> get gene clusters with Usearch
    

    # for command specific help
    MetaCHIP <command> -h
    '''

    print(help_message)


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help="--", dest='subparser_name')


    # arguments for PrepIn
    PI_parser = subparsers.add_parser('PI', description='Prepare input files', epilog='Example: MetaCHIP PI -h')
    PI_parser.add_argument('-i',                required=True, help='input genome folder')
    PI_parser.add_argument('-taxon',            required=True, help='taxonomic classification')
    PI_parser.add_argument('-p',                required=True, help='output prefix')
    PI_parser.add_argument('-r',                required=True, help='grouping rank')
    PI_parser.add_argument('-x',                required=False, default='fasta', help='file extension')
    PI_parser.add_argument('-grouping_only',    required=False, action="store_true", help='run grouping only, deactivate gene calling and phylogenetic tree building')
    PI_parser.add_argument('-nonmeta',          required=False, action="store_true", help='annotation as non-metagenome assembled genomes (non-MAG)')
    PI_parser.add_argument('-noblast',          required=False, action="store_true", help='not run all-vs-all blastn')
    PI_parser.add_argument('-t',                required=False, type=int, default=1, help='number of threads')
    PI_parser.add_argument('-force',            required=False, action="store_true", help='overwrite previous results')
    PI_parser.add_argument('-quiet',            required=False, action="store_true", help='not report progress')


    # arguments for BM approach
    BM_parser = subparsers.add_parser('BM', description='Best-match approach', epilog='Example: MetaCHIP BM -h')
    BM_parser.add_argument('-p',         required=True, help='output prefix')
    BM_parser.add_argument('-r',         required=True, help='grouping rank')
    BM_parser.add_argument('-g',         required=False, default=None, help='grouping file')
    BM_parser.add_argument('-cov',       required=False, type=int, default=70, help='coverage cutoff, deafult: 70')
    BM_parser.add_argument('-al',        required=False, type=int, default=200, help='alignment length cutoff, deafult: 200')
    BM_parser.add_argument('-flk',       required=False, type=int, default=10, help='the length of flanking sequences to plot (Kbp), deafult: 10')
    BM_parser.add_argument('-ip',        required=False, type=int, default=90, help='identity percentile cutoff, deafult: 90')
    BM_parser.add_argument('-ei',        required=False, type=float, default=95, help='end match identity cutoff, deafult: 95')
    BM_parser.add_argument('-t',         required=False, type=int, default=1, help='number of threads, deafult: 1')
    BM_parser.add_argument('-plot_iden', required=False, action="store_true", help='plot identity distribution')
    BM_parser.add_argument('-NoEbCheck', required=False, action="store_true", help='disable end break and contig match check')
    BM_parser.add_argument('-force',     required=False, action="store_true", help='overwrite previous results')
    BM_parser.add_argument('-quiet',     required=False, action="store_true", help='Do not report progress')
    BM_parser.add_argument('-tmp',       required=False, action="store_true", help='keep temporary files')


    # arguments for PG approach
    PG_parser = subparsers.add_parser('PG', description='Phylogenetic approach', epilog='Example: MetaCHIP PG -h')
    PG_parser.add_argument('-p',     required=True,  help='output prefix')
    PG_parser.add_argument('-r',     required=True,  help='grouping rank')
    PG_parser.add_argument('-g',     required=False, help='grouping file')
    PG_parser.add_argument('-tree',  required=False, default=None, help='species (SCG) tree')
    PG_parser.add_argument('-cov',   required=False, type=int, default=70, help='coverage cutoff, deafult: 70')
    PG_parser.add_argument('-al',    required=False, type=int, default=200, help='alignment length cutoff, deafult: 200')
    PG_parser.add_argument('-flk',   required=False, type=int, default=10, help='the length of flanking sequences to plot (Kbp), deafult: 10')
    PG_parser.add_argument('-ip',    required=False, type=int, default=90, help='identity percentile, deafult: 90')
    PG_parser.add_argument('-ei',    required=False, type=float, default=95, help='end match identity cutoff, deafult: 95')
    PG_parser.add_argument('-t',     required=False, type=int, default=1, help='number of threads, deafult: 1')
    PG_parser.add_argument('-force', required=False, action="store_true", help='overwrite previous results')
    PG_parser.add_argument('-quiet', required=False, action="store_true", help='Do not report progress')


    # get and check options
    args = None
    if (len(sys.argv) == 1) or (sys.argv[1] == '-h') or (sys.argv == '--help'):
        print_main_help()
        sys.exit(0)
    else:
        args = vars(parser.parse_args())


    if args['subparser_name'] == 'PI':
        PI(args)
    if args['subparser_name'] == 'BM':
        BM(args)
    if args['subparser_name'] == 'PG':
        PG(args)
