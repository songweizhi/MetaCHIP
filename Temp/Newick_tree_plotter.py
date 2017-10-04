from ete3 import Tree
from ete3 import TreeStyle
from ete3 import NodeStyle
from ete3 import TextFace


def plot_tree(tree_newick, tree_title):

    tree = Tree(tree_newick, format = 1)
    # set tree parameters
    ts = TreeStyle()
    ts.mode = "r"  # tree model: 'r' for rectangular, 'c' for circular
    ts.show_leaf_name = 0
    # set tree title text parameters
    ts.title.add_face(TextFace(tree_title,
                               fsize = 8,
                               fgcolor = 'black',
                               ftype = 'Arial',
                               tight_text = False),
                      column = 0)  # tree title text setting
    # set layout parameters
    ts.rotation = 0  # from 0 to 360
    ts.show_scale = False
    ts.margin_top = 10  # top tree image margin
    ts.margin_bottom = 10  # bottom tree image margin
    ts.margin_left = 10  # left tree image margin
    ts.margin_right = 10  # right tree image margin
    ts.show_border = False  # set tree image border
    ts.branch_vertical_margin = 3  # 3 pixels between adjancent branches

    # set tree node style
    for each_node in tree.traverse():
        # leaf node parameters
        if each_node.is_leaf():
            ns = NodeStyle()
            ns["shape"] = "circle"  # dot shape: circle, square or sphere
            ns["size"] = 0  # dot size
            ns['hz_line_width'] = 0.5  # branch line width
            ns['vt_line_width'] = 0.5  # branch line width
            ns['hz_line_type'] = 0  # branch line type: 0 for solid, 1 for dashed, 2 for dotted
            ns['vt_line_type'] = 0  # branch line type
            ns["fgcolor"] = "blue"  # the dot setting
            each_node.add_face(TextFace(each_node.name,
                                        fsize = 5,
                                        fgcolor = 'black',
                                        tight_text = False,
                                        bold = False),
                               column = 0,
                               position = 'branch-right')  # leaf node the node name text setting

            each_node.set_style(ns)

        # non-leaf node parameters
        else:
            nlns = NodeStyle()
            nlns["size"] = 0  # dot size
            #nlns["rotation"] = 45
            each_node.add_face(TextFace(each_node.name,

                                        fsize = 3,
                                        fgcolor = 'black',
                                        tight_text = False,
                                        bold = False),
                               column = 5,
                               position = 'branch-top')  # non-leaf node name text setting)

            each_node.set_style(nlns)

    tree.render('/Users/weizhisong/Desktop/Tree_' + tree_title + '.png', w = 900, units = "px", tree_style = ts)  # set figures size


tree_1 = '((((((((D5_amphiroa_bin16,D6_caulerpa_bin21)Saprospirace.)Chitinophaga.)Chitinophagi.,(((C7_amphiroa_bin21)Flavobacteri.)unclassified.,(((C1_amphiroa_bin24)Fluviicola)Crocinitomic.)Flavobacteri.)Flavobacteri.,((D18_ulva_bin9)Sphingobacte.)Sphingobacte.)Bacteroidete.)Bacteroidete.)FCB group,((((F9_ulva_bin45,F10_ulva_bin74)Rhodobactera.,((G5_delisea_bin38,G16_ulva_bin8,G7_delisea_bin9,G8_ecklonia_bin2,G10_ecklonia_bin7)Robiginitoma.)Hyphomonadac.)Rhodobactera.,G24_ulva_bin64,G25_ulva_bin65,G26_ulva_bin67,G18_ecklonia_bin26,G21_halophila_bin34,G23_halophila_bin48)Alphaproteob.,(L12_water_bin80,L6_ecklonia_bin49)Gammaproteob.)Proteobacter.,((((B2_delisea_bin83,B1_delisea_bin26)Planctomycet.)Planctomycet.)Planctomycet.,(((((A1_ecklonia_bin31)Rubritalea)Rubritaleace.)Verrucomicro.)Verrucomicro.)Verrucomicro.)PVC group)Bacteria);'


tree_2 = '(Refined_9:0.29446,Refined_34:0.25674,Refined_53:0.70393);'


plot_tree(tree_2, 'Katana')

