import os
from Bio import SeqIO
from Bio.Graphics import GenomeDiagram
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord


def plot_gbk(pwd_gbk_file, pwd_plot_file):
    contig = SeqIO.read(pwd_gbk_file, "genbank")
    diagram = GenomeDiagram.Diagram()

    # add a graph track...
    gene_content_track = diagram.new_track(1,
                                           greytrack = True,
                                           name = contig.name,
                                           greytrack_labels = True,
                                           greytrack_font = 'Helvetica',
                                           greytrack_fontsize = 12,
                                           height = 0.35,
                                           start = 0,
                                           end = len(contig),
                                           scale = True,
                                           scale_fontsize = 6,
                                           scale_ticks = 1,
                                           scale_smalltick_interval = 10000,
                                           scale_largetick_interval = 10000)

    feature_set = gene_content_track.new_set(type = 'feature')

    # add features to feature set
    for feature in contig.features:
        if feature.type == "CDS":
            feature_set.add_feature(feature,
                                    color = colors.lightblue,
                                    label = True,
                                    sigil = 'ARROW',
                                    arrowshaft_height = 0.7,
                                    arrowhead_length = 0.4,
                                    label_size = 10,
                                    label_position = "middle")

    # Finally draw it in both formats,
    diagram.draw(format = 'linear',
                 orientation = 'landscape',
                 tracklines = 0,
                 pagesize=(50 * cm, 10 * cm),
                 fragments = 1)

    diagram.write(pwd_plot_file, 'EPS')


def get_flanking_region(input_gbk_file, HGT_candidate, flanking_length):

    wd, gbk_file = os.path.split(input_gbk_file)
    new_gbk_file = '%s/%s_%sbp_temp.gbk' % (wd, HGT_candidate, flanking_length)
    new_gbk_final_file = '%s/%s_%sbp.gbk' % (wd, HGT_candidate, flanking_length)
    new_fasta_final_file = '%s/%s_%sbp.fasta' % (wd, HGT_candidate, flanking_length)
    output_plot = '%s/%s_%sbp.eps' % (wd, HGT_candidate, flanking_length)


    # get flanking range of candidate
    input_gbk = SeqIO.parse(input_gbk_file, "genbank")
    new_start = 0
    new_end = 0
    contig_length = 0
    for record in input_gbk:
        for gene in record.features:
            # get contig length
            if gene.type == 'source':
                contig_length = int(gene.location.end)
            # get new start and end points
            elif 'locus_tag' in gene.qualifiers:
                if gene.qualifiers['locus_tag'][0] == HGT_candidate:
                    # get new start point
                    new_start = gene.location.start - flanking_length
                    if new_start < 0:
                        new_start = 0
                    # get new end point
                    new_end = gene.location.end + flanking_length
                    if new_end > contig_length:
                        new_end = contig_length


    # get genes within flanking region
    keep_gene_list = []
    input_gbk = SeqIO.parse(input_gbk_file, "genbank")
    for record in input_gbk:
        for gene in record.features:
            if 'locus_tag' in gene.qualifiers:
                if (gene.location.start < new_start) and (gene.location.end >= new_start):
                    keep_gene_list.append(gene.qualifiers['locus_tag'][0])
                    new_start = gene.location.start
                elif (gene.location.start > new_start) and (gene.location.end < new_end):
                    keep_gene_list.append(gene.qualifiers['locus_tag'][0])
                elif (gene.location.start <= new_end) and (gene.location.end > new_end):
                    keep_gene_list.append(gene.qualifiers['locus_tag'][0])
                    new_end = gene.location.end


    # remove genes not in flanking region from gbk file
    input_gbk = SeqIO.parse(input_gbk_file, "genbank")
    new_gbk = open(new_gbk_file, 'w')
    for record in input_gbk:
        new_record_features = []
        for gene in record.features:
            if gene.type == 'source':
                new_record_features.append(gene)
            elif 'locus_tag' in gene.qualifiers:
                if gene.qualifiers['locus_tag'][0] in keep_gene_list:
                    #print(record.seq[gene.location.start:gene.location.end])
                    new_record_features.append(gene)
        record.features = new_record_features
        SeqIO.write(record, new_gbk, 'genbank')
    new_gbk.close()


    # remove sequences not in flanking region
    new_gbk_full_length = SeqIO.parse(new_gbk_file, "genbank")
    new_gbk_final = open(new_gbk_final_file, 'w')
    new_fasta_final = open(new_fasta_final_file, 'w')
    for record in new_gbk_full_length:
        # get new sequence
        new_seq = record.seq[new_start:new_end]
        new_contig_length = len(new_seq)
        new_record = SeqRecord(new_seq,
                               id = record.id,
                               name = record.name,
                               description = record.description,
                               annotations = record.annotations)

        # get new location
        new_record_features_2 = []
        for gene in record.features:
            if gene.type == 'source':
                gene_location_new = ''
                if gene.location.strand == 1:
                    gene_location_new = FeatureLocation(0, new_contig_length, strand=+1)
                if gene.location.strand == -1:
                    gene_location_new = FeatureLocation(0, new_contig_length, strand=-1)
                gene.location = gene_location_new
                new_record_features_2.append(gene)
            elif 'locus_tag' in gene.qualifiers:
                #print(record.seq[gene.location.start:gene.location.end])
                gene_location_new = ''
                if gene.location.strand == 1:
                    if gene.location.start - new_start < 0:
                        gene_location_new = FeatureLocation(0,
                                                            gene.location.end - new_start,
                                                            strand=+1)
                    else:
                        gene_location_new = FeatureLocation(gene.location.start - new_start,
                                                            gene.location.end - new_start,
                                                            strand=+1)
                if gene.location.strand == -1:
                    if gene.location.start - new_start < 0:
                        gene_location_new = FeatureLocation(0,
                                                            gene.location.end - new_start,
                                                            strand=-1)
                    else:
                        gene_location_new = FeatureLocation(gene.location.start - new_start,
                                                            gene.location.end - new_start,
                                                            strand=-1)
                gene.location = gene_location_new
                new_record_features_2.append(gene)
        new_record.features = new_record_features_2
        SeqIO.write(new_record, new_gbk_final, 'genbank')
        SeqIO.write(new_record, new_fasta_final, 'fasta')


    new_gbk_final.close()
    new_fasta_final.close()

    os.system('rm %s' % new_gbk_file)

    # plot
    plot_gbk(new_gbk_final_file, output_plot)


input_gbk_file = '/Users/songweizhi/Desktop/AAM.gbk'
get_flanking_region(input_gbk_file, 'AAM_02055', 2000)
