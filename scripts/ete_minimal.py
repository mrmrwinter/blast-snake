from ete3 import Tree, faces, AttrFace, TreeStyle, NodeStyle

def layout(node):
    if node.is_leaf():
        N = AttrFace("name", fsize=20)
        faces.add_face_to_node(N, node, 0, position="aligned")

def get_example_tree():

    # Set dashed blue lines in all leaves
    # nst1 = NodeStyle()
    # nst1["bgcolor"] = "LightSteelBlue"
    # nst2 = NodeStyle()
    # nst2["bgcolor"] = "Moccasin"
    # nst3 = NodeStyle()
    # nst3["bgcolor"] = "DarkSeaGreen"
    # nst4 = NodeStyle()
    # nst4["bgcolor"] = "Khaki"


    t = Tree("test.nwk")
    # for n in t.traverse(): # this removes branchlengths, but why??
    #     n.dist = 0

    # would something lke t.get_leaves_by_name(name) work? startswith?
    # n1 = t.get_common_ancestor("MGorilla", "LGorilla")
    # n1.set_style(nst1)
    # n2 = t.get_common_ancestor("Chimp", "Bonobo")
    # n2.set_style(nst2)
    # n3 = t.get_common_ancestor("c1", "c2", "c3")
    # n3.set_style(nst3)
    # n4 = t.get_common_ancestor("b3", "b4")
    # n4.set_style(nst4)
    ts = TreeStyle()
    ts.layout_fn = layout
    ts.show_leaf_name = False

    # midpoint root a tree
    root_point = t.get_midpoint_outgroup()
    t.set_outgroup(root_point)

    # ts.mode = "c" # this makes it a circle tree
    # ts.root_opening_factor = 1
    return t, ts

if __name__ == "__main__":
    t, ts = get_example_tree()
    # t.render("node_background.png", w=400, tree_style=ts)
    t.show(tree_style=ts)
