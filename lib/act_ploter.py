import os
from sys import stdout
from Bio import SeqIO
from Bio.SeqUtils import GC
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink


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



def get_gbk_blast_act(candidates_file, gbk_file, name_to_group_number_dict, path_to_output_act_folder):
    matches = open(candidates_file)
    total = 0
    for match in matches:
        total += 1

    matches = open(candidates_file)
    n = 1
    for match in matches:
        match = match.strip()
        genes = sorted(match.split('\t'))
        folder_name = '___'.join(genes)
        stdout.write("\rProcessing %dth of %d HGT candidates: %s" % (n, total, folder_name))
        os.mkdir('%s/%s' % (path_to_output_act_folder, folder_name))

        # Extract gbk and faa files
        records = SeqIO.parse(gbk_file, 'genbank')
        gene_gc_dict = {}
        for record in records:
            for gene_f in record.features :
                if 'locus_tag' in gene_f.qualifiers:
                    for gene_1 in genes :
                        if gene_1 in gene_f.qualifiers["locus_tag"]:
                            start = gene_f.location.start.position
                            end = gene_f.location.end.position
                            orf = record.seq[start :end]
                            gene_gc = GC(orf)
                            gene_gc_dict[gene_1] = float("{0:.2f}".format(gene_gc))
                            SeqIO.write(record, '%s/%s/%s.gbk' % (path_to_output_act_folder, folder_name, gene_1), 'genbank')
                            SeqIO.write(record, '%s/%s/%s.fasta' % (path_to_output_act_folder, folder_name, gene_1), 'fasta')

        # Run Blast
        prefix_c = '%s/%s' % (path_to_output_act_folder, folder_name)
        query_c = '%s/%s.fasta' % (prefix_c, genes[0])
        subject_c = '%s/%s.fasta' % (prefix_c, genes[1])
        output_c = '%s/%s.txt' % (prefix_c, folder_name)
        parameters_c_n = ' -evalue 1e-5 -outfmt 6 -task blastn'
        command_blast = 'blastn -query %s -subject %s -out %s%s' % (query_c, subject_c, output_c, parameters_c_n)
        os.system(command_blast)

        # read in gbk files
        matche_pair_list = []
        for each_gene in genes:
            path_to_gbk_file = '%s/%s/%s.gbk' % (path_to_output_act_folder, folder_name, each_gene)
            gene_contig = SeqIO.read(path_to_gbk_file, "genbank")
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
                                                            name = gene1_contig.name,
                                                            greytrack = True,
                                                            greytrack_labels = 1,
                                                            greytrack_font = 'Helvetica',
                                                            greytrack_fontsize = 12,
                                                            height = 0.35,
                                                            start = 0,
                                                            end = len(gene1_contig),
                                                            scale = True,
                                                            scale_fontsize = 6,
                                                            scale_ticks = 1,
                                                            scale_smalltick_interval = 10000,
                                                            scale_largetick_interval = 10000)
            # add gene content track for gene2_contig
            contig_2_gene_content_track = diagram.new_track(1,
                                                            name = gene2_contig.name,
                                                            greytrack = True,
                                                            greytrack_labels = 1,
                                                            greytrack_font = 'Helvetica',
                                                            greytrack_fontsize = 12,
                                                            height = 0.35,
                                                            start = 0,
                                                            end = len(gene2_contig),
                                                            scale = True,
                                                            scale_fontsize = 6,
                                                            scale_ticks = 1,
                                                            scale_smalltick_interval = 10000,
                                                            scale_largetick_interval = 10000)

            # add blank feature/graph sets to each track
            feature_sets_1 = contig_1_gene_content_track.new_set(type = 'feature')
            feature_sets_2 = contig_2_gene_content_track.new_set(type = 'feature')

            # add gene features to 2 blank feature sets
            set_contig_track_features(gene1_contig, name_to_group_number_dict, genes, feature_sets_1)
            set_contig_track_features(gene2_contig, name_to_group_number_dict, genes, feature_sets_2)

            ####################################### add crosslink from blast results #######################################

            path_to_blast_result = '%s/%s/%s.txt' % (path_to_output_act_folder, folder_name, folder_name)
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

            ############################################### Draw and Export ################################################

            diagram.draw(format = "linear",
                         orientation = "landscape",
                         pagesize = (75 * cm, 25 * cm),
                         fragments = 1,
                         start = 0,
                         end = max_len)

            diagram.write('%s/%s.eps' % (path_to_output_act_folder, folder_name), "EPS")
            os.system('rm -r %s/%s' % (path_to_output_act_folder, folder_name))

        n += 1
