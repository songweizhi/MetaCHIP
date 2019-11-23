import os
import argparse
from MetaCHIP.MetaCHIP_config import config_dict


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


circos_HGT_usage ='''
===================== circos_HGT example commands ====================

MetaCHIP circos_HGT -in SpongeEMP_p31_cir_plot_matrix.csv

======================================================================
'''


def circos_HGT(args, config_dict):

    matrix_in =         args['in']
    circos_HGT_R =      config_dict['circos_HGT_R']

    matrix_in_path, matrix_in_basename, matrix_in_extension = sep_path_basename_ext(matrix_in)
    pwd_outplot = '%s/%s.png' % (matrix_in_path, matrix_in_basename)

    os.system('Rscript %s -m %s -p %s' % (circos_HGT_R, matrix_in, pwd_outplot))


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='', add_help=True)

    parser.add_argument('-in', required=True,  help='input matrix')

    args = vars(parser.parse_args())

    circos_HGT(args, config_dict)
