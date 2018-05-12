#!/usr/bin/env python
import os
from Bio import SeqIO
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink
from reportlab.lib import colors
from reportlab.lib.units import cm


def set_contig_track_features(gene_contig, name_group_dict, candidate_list, feature_set):
    # add features to feature set
    for feature in gene_contig.features:
        if feature.type == "CDS":
            # define label color
            if feature.qualifiers['locus_tag'][0] in candidate_list:
                label_color = colors.blue
            else:
                label_color = colors.black

            # change gene name
            bin_name_gbk_split = feature.qualifiers['locus_tag'][0].split('_')
            bin_name_gbk = '_'.join(bin_name_gbk_split[:-1])
            feature.qualifiers['locus_tag'][0] = name_group_dict[bin_name_gbk] + '_' + bin_name_gbk_split[-1]

            # strands
            color = None
            label_angle = 0
            if feature.location.strand == 1:
                label_angle = 45
                color = colors.lightblue
            elif feature.location.strand == -1:
                label_angle = -225
                color = colors.lightgreen
            # add feature
            feature_set.add_feature(feature,
                                    color = color,
                                    label = True,
                                    sigil = 'ARROW',
                                    arrowshaft_height = 0.7,
                                    arrowhead_length = 0.4,
                                    label_color = label_color,
                                    label_size = 10,
                                    label_angle = label_angle,
                                    label_position = "middle")


file_1 = '16s.fasta'
file_2 = 'ctg.fasta'
path_to_blast_result = 'blast_results.txt'
plot_file_name = 'test'

# Run Blast
parameters_c_n = '-evalue 1e-5 -outfmt 6 -task blastn'
command_blast = 'blastn -query %s -subject %s -out %s %s' % (file_1, file_2, path_to_blast_result, parameters_c_n)
os.system(command_blast)

# read in fasta files
matche_pair_list = []
for each_file in [file_1, file_2]:
    gene_contig = SeqIO.read(each_file, "fasta")
    matche_pair_list.append(gene_contig)
bin_record_list = []
bin_record_list.append(matche_pair_list)

# create an empty diagram
diagram = GenomeDiagram.Diagram()

# add tracks to diagram
max_len = 0
for gene1_contig, gene2_contig in bin_record_list:
    # set diagram track length
    max_len = max(max_len, len(gene1_contig), len(gene2_contig))

    # add gene content track for gene1_contig
    contig_1_gene_content_track = diagram.new_track(1,
                                                    name=gene1_contig.name,
                                                    greytrack=True,
                                                    greytrack_labels=1,
                                                    greytrack_font='Helvetica',
                                                    greytrack_fontsize=12,
                                                    height=0.35,
                                                    start=0,
                                                    end=len(gene1_contig),
                                                    scale=True,
                                                    scale_fontsize=6,
                                                    scale_ticks=1,
                                                    scale_smalltick_interval=10000,
                                                    scale_largetick_interval=10000)

    # add gene content track for gene2_contig
    contig_2_gene_content_track = diagram.new_track(1,
                                                    name=gene2_contig.name,
                                                    greytrack=True,
                                                    greytrack_labels=1,
                                                    greytrack_font='Helvetica',
                                                    greytrack_fontsize=12,
                                                    height=0.35,
                                                    start=0,
                                                    end=len(gene2_contig),
                                                    scale=True,
                                                    scale_fontsize=6,
                                                    scale_ticks=1,
                                                    scale_smalltick_interval=10000,
                                                    scale_largetick_interval=10000)

    # add crosslink from blast results
    blast_results = open(path_to_blast_result)
    # parse blast results
    for each_line in blast_results:
        each_line_split = each_line.split('\t')
        query = each_line_split[0]
        target = each_line_split[1]
        identity = float(each_line_split[2])
        query_start = int(each_line_split[6])
        query_end = int(each_line_split[7])
        target_start = int(each_line_split[8])
        target_end = int(each_line_split[9])

        # use color to reflect identity
        color = colors.linearlyInterpolatedColor(colors.white, colors.red, 50, 100, identity)

        # determine which is which (query/target to contig_1/contig_2)
        # if query is contig_1
        if query == gene1_contig.name :
            link = CrossLink((contig_1_gene_content_track, query_start, query_end),
                             (contig_2_gene_content_track, target_start, target_end),
                             color = color,
                             border = color,
                             flip = False)
            diagram.cross_track_links.append(link)

        # if query is contig_2
        elif query == gene2_contig.name:
            link = CrossLink((contig_2_gene_content_track, query_start, query_end),
                             (contig_1_gene_content_track, target_start, target_end),
                             color = color,
                             border = color,
                             flip = False)
            diagram.cross_track_links.append(link)

    # Draw and Export
    diagram.draw(format="linear",
                 orientation="landscape",
                 pagesize=(75 * cm, 25 * cm),
                 fragments=1,
                 start=0,
                 end=max_len)

    diagram.write('%s.eps' % plot_file_name, "EPS")
