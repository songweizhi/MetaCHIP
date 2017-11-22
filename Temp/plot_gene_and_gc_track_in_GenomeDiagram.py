from Bio import SeqIO
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import test_GenomeDiagram, Diagram
from reportlab.lib import colors
from reportlab.lib.units import cm
# test_GenomeDiagram.py (https://github.com/biopython/biopython/blob/master/Tests/test_GenomeDiagram.py)
# download test_GenomeDiagram.py and copy it to folder "Bio/Graphics/GenomeDiagram/".
# /Library/Frameworks/Python.framework/Versions/3.4/lib/python3.4/site-packages/Bio/Graphics/GenomeDiagram/

contig = SeqIO.read('/Users/songweizhi/Desktop/test.gbk', "genbank")







diagram = GenomeDiagram.Diagram()

# add a graph track...
gc_content_track = diagram.new_track(1,
                                     greytrack = True,
                                     name = contig.name + " GC content",
                                     greytrack_labels = True,
                                     greytrack_font = 'Helvetica',
                                     greytrack_fontsize = 12,
                                     height = 0.35,
                                     start = 0,
                                     end = len(contig),
                                     scale = True,
                                     scale_fontsize = 12,
                                     scale_ticks = 0)

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


gc_content_set = gc_content_track.new_set(type = "graph")
feature_set = gene_content_track.new_set(type = 'feature')

step = len(contig)//500
gc_content_set.new_graph(test_GenomeDiagram.apply_to_window(contig.seq,
                                                            step,
                                                            test_GenomeDiagram.calc_gc_content, step),
                         'GC content',
                         style= 'bar',
                         color=colors.lightgreen,
                         altcolor=colors.darkseagreen)


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
             pagesize=(50 * cm, 20 * cm),
             fragments = 1)

diagram.write('/Users/songweizhi/Desktop/test_1.eps', 'EPS')
