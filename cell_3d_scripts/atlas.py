import csv
import math
from collections import deque
from collections.abc import Generator
from pathlib import Path
from typing import Optional

from brainglobe_utils.cells.cells import Cell

OTHER_OFFSET = 20000
OTHER_SUFFIX = "-other"
OUTSIDE_REGION_ID = 0


class AtlasNode:

    region_id: int = -1

    acronym: str = ""

    name: str = ""

    volume_voxels: int = 0

    volume_mm3: float = 0.0

    cell_count: int = 0

    parent: Optional["AtlasNode"] = None

    children: list["AtlasNode"] = None

    child_leaf: Optional["AtlasNode"] = None

    def __init__(self, region_id: int, acronym: str, name: str):
        self.region_id = region_id
        self.acronym = acronym
        self.name = name
        self.children = []

    @property
    def cell_density_voxels(self):
        return self.cell_count / self.volume_voxels

    @property
    def cell_density_mm3(self):
        return self.cell_count / self.volume_mm3

    def __str__(self):
        return (
            f"<Region {self.region_id} ({self.acronym}). Volume: {self.volume_voxels} voxels, "
            f"{self.volume_mm3:.5g} mm3. AtlasNode @ {id(self)}>"
        )

    def __repr__(self):
        return self.__str__()

    def post_order(self) -> Generator["AtlasNode"]:
        for child in self.children:
            yield from child.post_order()
        yield self


class AtlasTree:

    root: AtlasNode

    nodes_map: dict[int, AtlasNode]

    @classmethod
    def parse_vaa3d(cls, filename: str | Path) -> "AtlasTree":
        with open(filename) as fh:
            reader = csv.reader(fh, delimiter=",")
            lines = [[c.strip() for c in row] for row in reader]

        header = lines[0]

        id_col = header.index("ID")
        acronym_col = header.index("Acronym")
        tree_col = header.index("Structures")
        leaf_col = header.index("is_leaf")

        row = lines[1]
        root = AtlasNode(region_id=int(row[id_col]), acronym=row[acronym_col], name=row[tree_col])
        nodes = {root.region_id: root}
        is_leaves = {root.region_id: row[leaf_col].lower() == "leaf"}

        stack = deque([(0, root)])
        for row in lines[2:]:
            current_indent = 0
            for current_indent, val in enumerate(row[tree_col:]):
                if val:
                    break

            while stack[-1][0] >= current_indent:
                stack.pop()
            parent_indent, parent_node = stack[-1]

            if row[id_col]:
                node_id = int(row[id_col])
                acronym = row[acronym_col]
                name = row[tree_col + current_indent]
            else:
                node_id = parent_node.region_id + OTHER_OFFSET
                acronym = parent_node.acronym + OTHER_SUFFIX
                name = parent_node.name

            node = AtlasNode(region_id=node_id, acronym=acronym, name=name)
            parent_node.children.append(node)
            node.parent = parent_node

            if node.region_id in nodes:
                raise ValueError(f"Region {node.region_id} found twice")
            nodes[node.region_id] = node

            is_leaves[node.region_id] = row[leaf_col].lower() == "leaf"

            if not row[id_col]:
                parent_node.child_leaf = node

            stack.append((current_indent, node))

        tree = AtlasTree()
        tree.root = root
        tree.nodes_map = nodes

        for node in nodes.values():
            if bool(node.children) == is_leaves[node.region_id]:
                raise ValueError(f"Region {node.region_id} child nodes doesn't match its leaf designation")

        return tree

    def get_node_from_canonical_id(self, region_id: int) -> AtlasNode:
        nodes = self.nodes_map
        if region_id not in nodes:
            raise ValueError(f"Region {region_id} not found")

        node = nodes[region_id]
        if node.child_leaf is None:
            return node
        return node.child_leaf

    def read_regions_volumes(self, filename: str | Path) -> None:
        nodes = self.nodes_map

        with open(filename) as fh:
            reader = csv.reader(fh, delimiter=",")
            lines = [[c.strip() for c in row] for row in reader]

        header = lines[0]
        id_col = header.index("region_id")
        count_col = header.index("aggregated_voxel_count")
        mm3_col = header.index("aggregated_volume_mm3")

        for row in lines[1:]:
            node = nodes[int(row[id_col])]
            node.volume_voxels = int(row[count_col])
            node.volume_mm3 = float(row[mm3_col])

        for node in nodes.values():
            leaf = node.child_leaf
            if leaf is not None:
                normal_nodes = [n for n in node.children if n is not leaf]
                voxels = sum(n.volume_voxels for n in normal_nodes)
                mm3 = sum(n.volume_mm3 for n in normal_nodes)

                leaf.volume_voxels = node.volume_voxels - voxels
                leaf.volume_mm3 = node.volume_mm3 - mm3

        for node in nodes.values():
            child_voxels = sum(n.volume_voxels for n in node.children)
            child_mm3 = sum(n.volume_mm3 for n in node.children)

            if node.children and child_voxels != node.volume_voxels:
                raise ValueError(
                    f"Node's children's voxels sum of {child_voxels} doesn't match node's voxels for {node}"
                )
            if node.children and not math.isclose(child_mm3, node.volume_mm3):
                raise ValueError(
                    f"Node's children's volume (mm3) sum of {child_mm3} doesn't match node's mm3 volume for {node}"
                )

    def count_cells(self, cells: list[Cell], metadata_key: str = "region_id") -> int:
        nodes = self.nodes_map
        outside_cells = 0

        for cell in cells:
            region_id = cell.metadata[metadata_key]
            if region_id == OUTSIDE_REGION_ID:
                outside_cells += 1
            else:
                node = nodes[region_id]
                node.cell_count += 1

        for node in nodes.values():
            if node.cell_count and node.children:
                raise ValueError(f"Got cells for non-leaf node {node}")

        for node in self.root.post_order():
            node.cell_count = sum(n.cell_count for n in node.children)

        assert self.root.cell_count == len(cells) - outside_cells
        return outside_cells
