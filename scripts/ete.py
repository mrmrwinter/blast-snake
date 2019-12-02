from ete3 import Tree, TreeStyle, TextFace, faces, AttrFace, NodeStyle, PhyloTree




t = Tree(snakemake.input[0])

#### here to .....
#     ts = TreeStyle()
#     ts.mode = "c"
#     ts.show_leaf_name = True
#
#
#     ts.root_opening_factor = 1
#
#     ts.title.add_face(TextFace("If you can read this then TextFace works", fsize=20), column=0)
# get_farthest_oldest_node
### ...here, is not working as it should. Read the documentation and figure out how it works cause ete seems pretty good once it's properly jammed through.

t.render(snakemake.output[0])
