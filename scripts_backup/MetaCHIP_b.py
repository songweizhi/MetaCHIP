import sys
import argparse
from scripts.PrepIn import PrepIn
from scripts.BM import BM
from scripts.PG import PG

todo = '''


'''


def print_main_help():

    help_message = ''' 
        ...::: MetaCHIP :::...
        
    PrepIn  -> prepare input files 
    BM      -> Best-match approach 
    PG      -> Phylogenetic approach

    # for command specific help
    MetaCHIP <command> -h
    '''

    print(help_message)


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help="--", dest='subparser_name')


    # arguments for PrepIn
    PrepIn_parser = subparsers.add_parser('PrepIn', description='Prepare input files', epilog='Example: MetaCHIP PrepIn -h')
    PrepIn_parser.add_argument('-i', required=True, help='input genome folder')
    PrepIn_parser.add_argument('-taxon', required=True, help='taxonomic classification')
    PrepIn_parser.add_argument('-p', required=True, help='output prefix')
    PrepIn_parser.add_argument('-r', required=False, default='c', help='grouping rank')
    PrepIn_parser.add_argument('-x', required=False, default='fasta', help='file extension')
    PrepIn_parser.add_argument('-grouping_only', required=False, action="store_true", help='run grouping only, deactivate gene calling and phylogenetic tree building')
    PrepIn_parser.add_argument('-nonmeta', required=False, action="store_true", help='annotation as non-metagenome assembled genomes (non-MAG)')
    PrepIn_parser.add_argument('-noblast', required=False, action="store_true", help='not run all-vs-all blastn')
    PrepIn_parser.add_argument('-t', required=False, type=int, default=1, help='number of threads')
    PrepIn_parser.add_argument('-force', required=False, action="store_true", help='overwrite previous results')
    PrepIn_parser.add_argument('-quiet', required=False, action="store_true", help='not report progress')


    # arguments for BM approach
    BM_parser = subparsers.add_parser('BM', description='Best-match approach', epilog='Example: MetaCHIP BM -h')
    BM_parser.add_argument('-p', required=True, help='output prefix')
    BM_parser.add_argument('-g', required=False, default=None, help='grouping file')
    BM_parser.add_argument('-cov', required=False, type=int, default=70, help='coverage cutoff, deafult: 70')
    BM_parser.add_argument('-al', required=False, type=int, default=200, help='alignment length cutoff, deafult: 200')
    BM_parser.add_argument('-flk', required=False, type=int, default=10, help='the length of flanking sequences to plot (Kbp), deafult: 10')
    BM_parser.add_argument('-ip', required=False, type=int, default=90, help='identity percentile cutoff, deafult: 90')
    BM_parser.add_argument('-ei', required=False, type=float, default=95, help='end match identity cutoff, deafult: 95')
    BM_parser.add_argument('-t', required=False, type=int, default=1, help='number of threads, deafult: 1')
    BM_parser.add_argument('-plot_iden', required=False, action="store_true", help='plot identity distribution')
    BM_parser.add_argument('-NoEbCheck', required=False, action="store_true", help='disable end break and contig match check')
    BM_parser.add_argument('-force', required=False, action="store_true", help='overwrite previous results')
    BM_parser.add_argument('-quiet', required=False, action="store_true", help='Do not report progress')
    BM_parser.add_argument('-tmp', required=False, action="store_true", help='keep temporary files')


    # arguments for PG approach
    PG_square_parser = subparsers.add_parser('PG', description='Phylogenetic approach', epilog='Example: MetaCHIP PG -h')
    PG_square_parser.add_argument('-p', required=True, help='output prefix')
    PG_square_parser.add_argument('-g', required=False, help='grouping file')
    PG_square_parser.add_argument('-o', required=True, help='orthologs folder')
    PG_square_parser.add_argument('-tree', required=False, default=None, help='species (SCG) tree')
    PG_square_parser.add_argument('-cov', required=False, type=int, default=70, help='coverage cutoff, deafult: 70')
    PG_square_parser.add_argument('-al', required=False, type=int, default=200, help='alignment length cutoff, deafult: 200')
    PG_square_parser.add_argument('-flk', required=False, type=int, default=10, help='the length of flanking sequences to plot (Kbp), deafult: 10')
    PG_square_parser.add_argument('-ip', required=False, type=int, default=90, help='identity percentile, deafult: 90')
    PG_square_parser.add_argument('-ei', required=False, type=float, default=95, help='end match identity cutoff, deafult: 95')
    PG_square_parser.add_argument('-t', required=False, type=int, default=1, help='number of threads, deafult: 1')
    PG_square_parser.add_argument('-force', required=False, action="store_true", help='overwrite previous results')
    PG_square_parser.add_argument('-quiet', required=False, action="store_true", help='Do not report progress')


    # get and check options
    args = None
    if (len(sys.argv) == 1) or (sys.argv[1] == '-h') or (sys.argv == '--help'):
        print_main_help()
        sys.exit(0)
    else:
        args = vars(parser.parse_args())


    if args['subparser_name'] == 'PrepIn':
        PrepIn(args)
    if args['subparser_name'] == 'BM':
        BM(args)
    if args['subparser_name'] == 'PG':
        PG(args)
