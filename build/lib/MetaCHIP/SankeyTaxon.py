import os
import argparse


SankeyTaxon_parser_usage ='''
===================================== SankeyTaxon example commands ====================================

# Visualize taxonomic classification at phylum, class and genus level
MetaCHIP SankeyTaxon -taxon KelpBins_gtdbtk.tsv -r pcg -p KelpBins

# Visualize taxonomic classification at phylum and genus level, with unclear classifications excluded 
MetaCHIP SankeyTaxon -taxon KelpBins_gtdbtk.tsv -r pg -p KelpBins -ec

=======================================================================================================
'''


def unique_list_elements(list_input):

    list_output = []
    for each_element in list_input:
        if each_element not in list_output:
            list_output.append(each_element)

    return list_output


def SankeyTaxon(args):

    GTDB_output =                       args['taxon']
    taxon_ranks =                       args['r']
    output_prefix =                     args['p']
    plot_width =                        args['x']
    plot_height =                       args['y']
    explicit_classifications_only =     args['ec']


    pwd_self = os.path.realpath(__file__)
    self_path = '/'.join(pwd_self.split('/')[:-1])
    pwd_get_sankey_plot_R = '%s/get_sankey_plot.R' % self_path

    output_file_txt = '%s_sankey_taxon.txt' % output_prefix

    rank_full_list = ['d', 'p', 'c', 'o', 'f', 'g', 's']
    rank_pos_dict = {'d': 0, 'p': 1, 'c': 2, 'o': 3, 'f': 4, 'g': 5, 's': 6}

    ranks_to_plot = []
    for each_rank in rank_full_list:
        if each_rank in taxon_ranks:
            if each_rank not in ranks_to_plot:
                ranks_to_plot.append(each_rank)


    paired_taxon_list_all = []
    genome_number = 0
    for each_genome_taxon in open(GTDB_output):
        if not each_genome_taxon.startswith('user_genome'):

            genome_number += 1

            genome_taxon_split = each_genome_taxon.strip().split('\t')

            taxon_asign_full = []
            if len(genome_taxon_split) == 1:
                taxon_asign_full = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']

            elif (len(genome_taxon_split) > 1) and (';' in genome_taxon_split[1]):
                assignment = genome_taxon_split[1].split(';')
                if len(assignment) == 7:
                    taxon_asign_full = assignment
                if len(assignment) == 6:
                    taxon_asign_full = assignment + ['s__']
                if len(assignment) == 5:
                    taxon_asign_full = assignment + ['g__', 's__']
                if len(assignment) == 4:
                    taxon_asign_full = assignment + ['f__', 'g__', 's__']
                if len(assignment) == 3:
                    taxon_asign_full = assignment + ['o__', 'f__', 'g__', 's__']
                if len(assignment) == 2:
                    taxon_asign_full = assignment + ['c__', 'o__', 'f__', 'g__', 's__']

            elif (len(genome_taxon_split) > 1) and (';' not in genome_taxon_split[1]):
                taxon_asign_full = [genome_taxon_split[1]] + ['p__', 'c__', 'o__', 'f__', 'g__', 's__']


            rank_to_plot_paired = []
            n = 0
            while n < len(ranks_to_plot) -1:
                rank_to_plot_paired.append([ranks_to_plot[n], ranks_to_plot[n + 1]])
                n += 1


            for each_pair in rank_to_plot_paired:

                each_pair_taxon_left = taxon_asign_full[rank_pos_dict[each_pair[0]]]
                each_pair_taxon_right = taxon_asign_full[rank_pos_dict[each_pair[1]]]
                each_pair_taxon = '%s,%s' % (each_pair_taxon_left, each_pair_taxon_right)

                # combine unclear classifications
                if explicit_classifications_only is False:
                    paired_taxon_list_all.append(each_pair_taxon)

                # don't show unclear classification
                else:
                    if (each_pair_taxon_left.endswith('__') is False) and (each_pair_taxon_right.endswith('__') is False):
                        paired_taxon_list_all.append(each_pair_taxon)


    paired_taxon_list_uniq = unique_list_elements(paired_taxon_list_all)
    paired_taxon_list_uniq_count_dict = {}
    for each_key in paired_taxon_list_uniq:
        paired_taxon_list_uniq_count_dict[each_key] = paired_taxon_list_all.count(each_key)


    output_file_handle = open(output_file_txt, 'w')
    output_file_handle.write('C1,C2,Number\n')
    for each_genome_taxon in paired_taxon_list_uniq_count_dict:
        for_out = '%s,%s\n' % (each_genome_taxon, paired_taxon_list_uniq_count_dict[each_genome_taxon])
        output_file_handle.write(for_out)
    output_file_handle.close()


    # run Rscript

    if plot_width == None:
        if len(ranks_to_plot) <= 5:
            plot_width = (len(ranks_to_plot) - 1) * 300
        else:
            plot_width = (len(ranks_to_plot)-1)*250

    if plot_height == None:
        plot_height = genome_number * 2
        if plot_height < 600:
            plot_height = 600

    os.system('Rscript %s -f %s -x %s -y %s' % (pwd_get_sankey_plot_R, output_file_txt, plot_width, plot_height))


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='', add_help=True)

    parser.add_argument('-taxon', required=True,  help='taxon classification results')
    parser.add_argument('-r',     required=True,  help='taxon ranks to plot, e.g. dpcofgs, pco, pcf, cfs')
    parser.add_argument('-p',     required=True,  help='output prefix')
    parser.add_argument('-ec',    required=False, action='store_true', help='only plot explicit classifications')
    parser.add_argument('-x',     required=False, type=int, help='plot width')
    parser.add_argument('-y',     required=False, type=int, help='plot height')

    args = vars(parser.parse_args())

    SankeyTaxon(args)

