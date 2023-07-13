import pytest

import anytree
import anytree.exporter

from sprog import tree, utils


def test_load_and_save_tree():
    outprefix = "tmp.test_tree"
    utils.syscall(f"rm -rf {outprefix}*")
    tmp_newick = f"{outprefix}.newick"
    tmp_tsv = f"{outprefix}.tsv"
    tmp_json = f"{outprefix}.json"

    # Tree looks like this:
    #                ┌─subsp1
    #     ┌─speciesA─┤
    #     │          └─subsp2
    # ────┤
    #     ├─speciesB
    #     └─speciesC
    expect_tree_root = anytree.Node("root")
    speciesA = anytree.Node("speciesA", parent=expect_tree_root)
    anytree.Node("speciesB", parent=expect_tree_root)
    anytree.Node("speciesC", parent=expect_tree_root)
    anytree.Node("subsp1", parent=speciesA)
    anytree.Node("subsp2", parent=speciesA)
    exporter = anytree.exporter.JsonExporter(indent=2, sort_keys=True)
    expect_json = exporter.export(expect_tree_root)
    expect_nodes = ["root", "speciesA", "speciesB", "speciesC", "subsp1", "subsp2"]

    with open(tmp_newick, "w") as f:
        print("((subsp1,subsp2)speciesA,speciesB,speciesC);", file=f)

    with open(tmp_tsv, "w") as f:
        print("speciesA", "subsp1", sep="\t", file=f)
        print("speciesA", "subsp2", sep="\t", file=f)
        print("speciesB", sep="\t", file=f)
        print("speciesC", sep="\t", file=f)

    t = tree.Tree(tmp_newick)
    assert t.to_json(filename=tmp_json) == expect_json
    assert sorted(list(t.nodes.keys())) == expect_nodes
    assert [t.nodes[x].name for x in sorted(t.nodes)] == expect_nodes

    t = tree.Tree(tmp_tsv)
    assert t.to_json() == expect_json
    assert sorted(list(t.nodes.keys())) == expect_nodes
    assert [t.nodes[x].name for x in sorted(t.nodes)] == expect_nodes

    t = tree.Tree(tmp_json)
    assert t.to_json() == expect_json
    assert sorted(list(t.nodes.keys())) == expect_nodes
    assert [t.nodes[x].name for x in sorted(t.nodes)] == expect_nodes

    assert t.find_node_by_name("speciesC").name == "speciesC"
    assert t.find_node_by_name("not_in_tree") is None

    utils.syscall(f"rm -f {outprefix}*")


def test_leaf_combo_parent_nodes():
    outprefix = "tmp.test_leaf_combo_parent_nodes"
    utils.syscall(f"rm -rf {outprefix}*")
    tmp_newick = f"{outprefix}.newick"
    with open(tmp_newick, "w") as f:
        print(
            "(((subsp1,subsp2,subsp3)speciesA,speciesB,speciesC)genus1,(speciesX,speciesY)genus2);",
            file=f,
        )
    t = tree.Tree(tmp_newick)
    # print(anytree.RenderTree(t.tree_root, style=anytree.render.ContStyle()))
    # tree looks like this:
    #    Node('/root')
    # ├── Node('/root/genus1')
    # │   ├── Node('/root/genus1/speciesA')
    # │   │   ├── Node('/root/genus1/speciesA/subsp1')
    # │   │   ├── Node('/root/genus1/speciesA/subsp2')
    # │   │   └── Node('/root/genus1/speciesA/subsp3')
    # │   ├── Node('/root/genus1/speciesB')
    # │   └── Node('/root/genus1/speciesC')
    # └── Node('/root/genus2')
    #    ├── Node('/root/genus2/speciesX')
    #    └── Node('/root/genus2/speciesY')

    leaf_names = [
        "subsp1",
        "subsp2",
        "subsp3",
        "speciesB",
        "speciesC",
        "speciesX",
        "speciesY",
    ]
    t.init_leaf_combination_data(leaf_names)
    # check leaf combinations correct. The values are actually Node
    # objects. We check the names match instead of creating nodes
    expect = [
        ({0, 1, 2, 3, 4, 5, 6}, "root"),
        ({0, 1, 2, 3, 4}, "genus1"),
        ({0, 1, 2}, "speciesA"),
        ({5, 6}, "genus2"),
    ]
    assert sorted(expect) == sorted((x[0], x[1].name) for x in t.leaf_combos)

    got = t.leaf_combo_parent_nodes({0}, min_prop_contain=0.95, max_prop_outside=0.05)
    assert got.name == "subsp1"
    got = t.leaf_combo_parent_nodes({1}, min_prop_contain=0.95, max_prop_outside=0.05)
    assert got.name == "subsp2"
    got = t.leaf_combo_parent_nodes({2}, min_prop_contain=0.95, max_prop_outside=0.05)
    assert got.name == "subsp3"
    got = t.leaf_combo_parent_nodes(
        {0, 1}, min_prop_contain=0.95, max_prop_outside=0.05
    )
    assert got is None
    got = t.leaf_combo_parent_nodes(
        {0, 1}, min_prop_contain=0.68, max_prop_outside=0.05
    )
    assert got is None
    got = t.leaf_combo_parent_nodes(
        {0, 1}, min_prop_contain=0.66, max_prop_outside=0.05
    )
    assert got.name == "speciesA"
    got = t.leaf_combo_parent_nodes(
        {0, 1, 4}, min_prop_contain=0.66, max_prop_outside=0.05
    )
    assert got is None
    got = t.leaf_combo_parent_nodes(
        {0, 1, 4}, min_prop_contain=0.66, max_prop_outside=0.5
    )
    assert got.name == "speciesA"
    utils.syscall(f"rm -f {outprefix}*")


def test_leaf_combinations_to_nodes():
    outprefix = "tmp.test_leaf_combos"
    utils.syscall(f"rm -rf {outprefix}*")
    tmp_newick = f"{outprefix}.newick"
    with open(tmp_newick, "w") as f:
        print("((subsp1,subsp2)speciesA,speciesB,speciesC);", file=f)
    t = tree.Tree(tmp_newick)
    leaf_names = ["subsp1", "subsp2", "speciesB", "speciesC"]
    got = t.leaf_combinations_to_nodes(leaf_names)
    assert len(got) == 6
    assert got[(0,)].name == "subsp1"
    assert got[(1,)].name == "subsp2"
    assert got[(0, 1)].name == "speciesA"
    assert got[(2,)].name == "speciesB"
    assert got[(3,)].name == "speciesC"
    assert got[(0, 1, 2, 3)].name == "root"
    utils.syscall(f"rm -f {outprefix}*")
