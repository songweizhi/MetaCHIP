import os
from ete3 import Tree
from ete3 import TreeStyle
from ete3 import NodeStyle
from ete3 import TextFace


def plot_tree(tree, tree_title, tree_output):
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

    tree.render(tree_output, w=900, units="px", tree_style=ts)  # set figures size


os.chdir('/Users/songweizhi/Desktop')

tree = Tree('/Users/songweizhi/Desktop/species_tree.newick', format = 0)
plot_tree(tree, 'Overall', 'Overall.png')

sub_tree = Tree('/Users/songweizhi/Desktop/species_tree.newick', format = 0)
sub_tree.prune(['AMS', 'AKV', 'BDS', 'BAD'])
plot_tree(sub_tree, 'sub_tree', 'sub_tree.png')
