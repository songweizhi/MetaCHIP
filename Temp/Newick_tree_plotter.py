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


#os.chdir('/Users/songweizhi/Desktop')

tree_2 = '(CF_Refined_71:0.21847,CF_Refined_170:0.41504,(((CF_Refined_7:0.63495,CF_Refined_96:0.68718)0.984:0.33915,CF_Refined_82:0.16074)0.980:0.12012,((CF_Refined_25:0.20437,CF_Refined_95:1.40476)0.367:0.60450,(((CF_Refined_86:0.37933,(CF_Refined_74:0.61406,CF_Refined_43:0.10850)1.000:0.34468)0.999:0.19175,((CF_Refined_57:0.19003,(CF_Refined_99:0.18534,CF_Refined_160:0.33153)0.861:0.04660)1.000:0.18553,(CF_Refined_78:0.64747,(CF_Refined_129:0.26317,(CF_Refined_64:0.25293,CF_Refined_100:0.10449)0.993:0.14949)0.577:0.13006)1.000:0.19231)0.998:0.16057)0.961:0.08413,((CF_Refined_155:0.18262,CF_Refined_180:0.02711)1.000:0.42318,(CF_Refined_162:0.64092,CF_Refined_131:0.36385)1.000:0.48656)0.992:0.23471)0.853:0.05581)0.955:0.08666)1.000:0.14706);'

tree = Tree(tree_2, format = 1)


plot_tree(tree, 'Species_Tree', 'Species_Tree.png')

