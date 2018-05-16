from ete3 import Tree, TreeStyle, NodeStyle, TextFace


# tree: tree in newick format
# tree type: species tree or gene tree
# name_list: if one node in name_list, it's name will be displayed in red
def plot_gene_tree(tree, tree_type, gene_name, tree_file_name, name_list, tree_image_folder):
    # set tree parameters
    ts = TreeStyle()
    ts.mode = "r"  # tree model: 'r' for rectangular, 'c' for circular
    #ts.scale = 50
    ts.show_leaf_name = False
    tree_title = '%s (%s)' % (tree_type, gene_name)  # define tree title
    # tree title text setting
    ts.title.add_face(TextFace(tree_title,
                               fsize = 8,
                               fgcolor = 'black',
                               ftype = 'Arial',
                               tight_text = False),
                      column = 0)

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
        if each_node.is_leaf():  # leaf node parameters
            ns = NodeStyle()
            ns["shape"] = "circle"  # dot shape: circle, square or sphere
            ns["size"] = 0  # dot size
            ns['hz_line_width'] = 0.5  # branch line width
            ns['vt_line_width'] = 0.5  # branch line width
            ns['hz_line_type'] = 0  # branch line type: 0 for solid, 1 for dashed, 2 for dotted
            ns['vt_line_type'] = 0  # branch line type
            if each_node.name in name_list:
                ns["fgcolor"] = "red"  # the dot setting
                # the node name text setting
                each_node.add_face(TextFace(each_node.name,
                                            fsize = 8,
                                            fgcolor = 'red',
                                            tight_text = False,
                                            bold = False),
                                   column = 0,
                                   position = 'branch-right')
                each_node.set_style(ns)
            else:
                ns["fgcolor"] = "blue"  # the dot setting
                # the node name text setting
                each_node.add_face(TextFace(each_node.name,
                                            fsize = 8,
                                            fgcolor = 'black',
                                            tight_text = False,
                                            bold = False),
                                   column = 0,
                                   position = 'branch-right')
                each_node.set_style(ns)
        else:  # non-leaf node parameters
            nlns = NodeStyle()
            nlns["size"] = 0  # dot size
            each_node.set_style(nlns)
    # set figures size
    tree.render('%s/%s_%s.png' % (tree_image_folder, tree_type, tree_file_name), w = 900, units = "px", tree_style = ts)


def plot_species_tree(tree_newick, tree_type, gene_name, tree_file_name, name_list, tree_image_folder):
    # set tree parameters
    tree = Tree(tree_newick, format = 8)
    ts = TreeStyle()
    ts.mode = "r"  # tree model: 'r' for rectangular, 'c' for circular
    ts.show_leaf_name = False
    tree_title = tree_type + ' (' + gene_name + ')'  # define tree title
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
            ns['shape'] = 'circle'  # dot shape: circle, square or sphere
            ns['size'] = 0  # dot size
            ns['hz_line_width'] = 0.5  # branch line width
            ns['vt_line_width'] = 0.5  # branch line width
            ns['hz_line_type'] = 0  # branch line type: 0 for solid, 1 for dashed, 2 for dotted
            ns['vt_line_type'] = 0  # branch line type
            if each_node.name in name_list:
                ns['fgcolor'] = 'red'  # the dot setting
                each_node.add_face(TextFace(each_node.name,
                                            fsize = 8,
                                            fgcolor = 'red',
                                            tight_text = False,
                                            bold = False),
                                   column = 0,
                                   position = 'branch-right')  # the node name text setting
                each_node.set_style(ns)
            else:
                ns['fgcolor'] = 'blue'  # the dot setting
                each_node.add_face(TextFace(each_node.name,
                                            fsize = 8,
                                            fgcolor = 'black',
                                            tight_text = False,
                                            bold = False),
                                   column = 0,
                                   position = 'branch-right')  # the node name text setting
                each_node.set_style(ns)

        # non-leaf node parameters
        else:
            nlns = NodeStyle()
            nlns['size'] = 0  # dot size
            each_node.add_face(TextFace(each_node.name,
                                        fsize = 4,
                                        fgcolor = 'black',
                                        tight_text = False,
                                        bold = False),
                               column = 5,
                               position = 'branch-top')  # non-leaf node name text setting)
            each_node.set_style(nlns)
    # set figures size
    tree.render('%s/%s_%s.png' % (tree_image_folder, tree_type, tree_file_name), w = 900, units = 'px', tree_style = ts)

