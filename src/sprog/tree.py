import logging

import anytree
import anytree.exporter
import anytree.importer
import newick

from sprog import utils

EDGES_TO_MYK_PANEL_TYPE = [
    "complex",
    "sub-complex",
    "species",
    "sub-species",
]


class Tree:
    def __init__(self, infile, file_format=None):
        self.tree_root = anytree.Node("root")
        self.nodes = {"root": self.tree_root}
        if file_format is None:
            file_format = self.guess_tree_file_format(infile)
        logging.debug(f"Guessed format of tree file '{infile}' as '{file_format}'")

        if file_format not in ["newick", "tsv", "json"]:
            raise Exception(f"Unrecognised file format: {file_format}")
        exec(f"self.load_from_{file_format}(infile)")

    def guess_tree_file_format(self, filename):
        with open(filename) as f:
            lines = [x.rstrip() for x in f.readlines()]

        if lines[0].strip().startswith("(") and lines[-1][-1] == ";":
            return "newick"
        elif lines[0].strip().startswith("{") and lines[-1][-1] == "}":
            return "json"
        else:
            return "tsv"

    def load_from_newick(self, filename):
        newick_tree = newick.read(filename)
        logging.debug(
            "loaded from newick file:\n" + str(newick_tree[0].ascii_art()) + "\n"
        )
        for node in newick_tree[0].walk():
            # if this is the root node, force its name to "root" and skip it
            if node.ancestor is None:
                node.name = "root"
                continue

            assert node.name is not None
            assert node.name not in self.nodes  # names should be unique
            # get ancestor and check it is already in the tree
            assert node.ancestor is not None
            assert node.ancestor.name in self.nodes
            new_node = anytree.Node(node.name, parent=self.nodes[node.ancestor.name])
            self.nodes[node.name] = new_node

    def load_from_tsv(self, filename):
        with open(filename) as f:
            for line in map(str.rstrip, f):
                fields = line.split("\t")
                parent = self.nodes["root"]
                for name in fields:
                    if name not in self.nodes:
                        self.nodes[name] = anytree.Node(name, parent=parent)
                    parent = self.nodes[name]

    def to_json(self, filename=None):
        exporter = anytree.exporter.JsonExporter(indent=2, sort_keys=True)
        json_str = exporter.export(self.tree_root)
        if filename is not None:
            with open(filename, "w") as f:
                print(json_str, file=f)
        return json_str

    def load_from_json(self, filename):
        importer = anytree.importer.JsonImporter()
        with open(filename) as f:
            self.tree_root = importer.read(f)
        self.nodes = {}
        for node in anytree.PreOrderIter(self.tree_root):
            assert node.name not in self.nodes
            self.nodes[node.name] = node

    def init_leaf_combination_data(self, leaf_names_list):
        self.leaf_combo_name2index = {name: i for i, name in enumerate(leaf_names_list)}
        self.leaf_combo_index2leaf = {i: name for i, name in enumerate(leaf_names_list)}
        self.leaf_combos = []
        used_leaves = set()

        for node in anytree.PreOrderIter(self.tree_root):
            if node.is_leaf:
                index = self.leaf_combo_name2index[node.name]
                self.leaf_combo_index2leaf[index] = node
                used_leaves.add(node.name)
                continue

            leaf_combo = set(
                self.leaf_combo_name2index[n.name] for n in node.leaves if n.is_leaf
            )
            if len(leaf_combo) == 1:
                continue

            self.leaf_combos.append((leaf_combo, node))

        assert set(leaf_names_list) == used_leaves

    def leaf_combo_parent_nodes(
        self, leaf_combo, min_prop_contain=0.95, max_prop_outside=0.05
    ):
        if len(leaf_combo) == 1:
            return self.leaf_combo_index2leaf[list(leaf_combo)[0]]

        for test_combo, combo_parent_node in self.leaf_combos:
            intersection = leaf_combo.intersection(test_combo)
            if len(intersection) < min_prop_contain * len(test_combo):
                continue

            non_test_combo_leaves = len(self.leaf_combo_name2index) - len(test_combo)
            if non_test_combo_leaves == 0:
                continue

            not_intersect = len(leaf_combo) - len(intersection)
            if not_intersect <= max_prop_outside * non_test_combo_leaves:
                return combo_parent_node

        return None

    def leaf_combinations_to_nodes(self, leaf_names):
        used_leaves = set()
        leaf_combos = {}

        for node in anytree.PreOrderIter(self.tree_root):
            leaf_nodes = set(n.name for n in node.leaves if n.is_leaf)
            used_leaves.update(leaf_nodes)
            combo = tuple(i for i, name in enumerate(leaf_names) if name in leaf_nodes)
            leaf_combos[combo] = node

        assert used_leaves == set(leaf_names)
        return leaf_combos

    def write_to_file(self, outfile, node_unitig_combos=None):
        if node_unitig_combos is None:
            node_unitig_combos = {}
        l = lambda n: n.name + " " + str(len(node_unitig_combos.get(n, [])))
        with open(outfile, "w") as f:
            print(
                anytree.RenderTree(self.tree_root, style=anytree.ContStyle()).by_attr(
                    l
                ),
                file=f,
            )

    def find_node_by_name(self, name):
        return anytree.search.find(self.tree_root, lambda n: n.name == name)

    def edges_to_root(self, node):
        if type(node) == str:
            node = self.find_node_by_name()
        return node.depth

    def myk_panel_type(self, node):
        return EDGES_TO_MYK_PANEL_TYPE[self.edges_to_root(node)]

    def write_mykrobe_hierarchy_json(self, outfile):
        harch = {
            "Non_tuberculosis_mycobacterium_complex": {
                "phylo_group": "complex",
                "children": {
                    "subNon_tuberculosis_mycobacterium_complex": {
                        "phylo_group": "sub-complex",
                        "children": {},
                    },
                },
            },
            "Mycobacterium_tuberculosis_complex": {
                "phylo_group": "complex",
                "children": {
                    "subMycobacterium_tuberculosis_complex": {
                        "phylo_group": "sub-complex",
                        "children": {},
                    },
                },
            },
        }

        complex_lookup = {
            "subNon_tuberculosis_mycobacterium_complex": harch[
                "Non_tuberculosis_mycobacterium_complex"
            ]["children"]["subNon_tuberculosis_mycobacterium_complex"],
            "subMycobacterium_tuberculosis_complex": harch[
                "Mycobacterium_tuberculosis_complex"
            ]["children"]["subMycobacterium_tuberculosis_complex"],
        }

        for node in anytree.PreOrderIter(self.tree_root):
            if node.depth == 0:
                pass
            elif node.depth == 1:
                assert node.name in complex_lookup
            elif node.depth == 2:
                parent_complex = complex_lookup[node.parent.name]["children"]
                name = "_".join(node.name.split())
                assert name not in parent_complex
                parent_complex[name] = {
                    "phylo_group": "species",
                    "children": {},
                }
            else:
                assert node.depth == 3
                species_node = node.parent
                complex_node = species_node.parent
                species_name = "_".join(species_node.name.split())
                species_entry = complex_lookup[complex_node.name]["children"][
                    species_name
                ]["children"]
                name = "_".join(node.name.split())
                assert name not in species_entry
                species_entry[name] = {
                    "phylo_group": "sub-species",
                    "children": {},
                }

        utils.write_json(harch, outfile)
