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
        if self.volume_voxels:
            return self.cell_count / self.volume_voxels
        return 0

    @property
    def cell_density_mm3(self):
        if self.volume_mm3:
            return self.cell_count / self.volume_mm3
        return 0

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

    def pre_order(self) -> Generator["AtlasNode"]:
        yield self
        for child in self.children:
            yield from child.pre_order()

    def pre_order_with_depth(self, depth=0) -> Generator[tuple["AtlasNode", int]]:
        yield self, depth
        for child in self.children:
            yield from child.pre_order_with_depth(depth + 1)


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
        acronym_col = header.index("Abbreviation")
        tree_col = header.index("RegionLevel_01")
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

            is_other_node = row[acronym_col].endswith(OTHER_SUFFIX)
            if is_other_node:
                node_id = parent_node.region_id + OTHER_OFFSET
                acronym = parent_node.acronym + OTHER_SUFFIX

                if node_id != int(row[id_col]):
                    raise ValueError(f"Expected the {OTHER_SUFFIX} region to be {node_id}, got {row[id_col]}")
                if acronym != row[acronym_col]:
                    raise ValueError(f"Expected the {OTHER_SUFFIX} region to be {acronym}, got {row[acronym_col]}")

            node_id = int(row[id_col])
            acronym = row[acronym_col]
            name = row[tree_col + current_indent]

            node = AtlasNode(region_id=node_id, acronym=acronym, name=name)
            parent_node.children.append(node)
            node.parent = parent_node

            if node.region_id in nodes:
                raise ValueError(f"Region {node.region_id} found twice")
            nodes[node.region_id] = node

            is_leaves[node.region_id] = row[leaf_col].lower() == "leaf"

            if is_other_node:
                parent_node.child_leaf = node

            stack.append((current_indent, node))

        tree = AtlasTree()
        tree.root = root
        tree.nodes_map = nodes

        for node in nodes.values():
            if bool(node.children) == is_leaves[node.region_id]:
                raise ValueError(f"Region {node.region_id} child nodes doesn't match its leaf designation")

        if len(nodes) != len(lines) - 1:
            raise ValueError(f"Expected {len(lines) - 1} unique tree nodes, got {len(nodes)}")

        return tree

    def get_leaf_node_from_canonical_id(self, region_id: int) -> AtlasNode:
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
            if node.children:
                node.cell_count = sum(n.cell_count for n in node.children)

        assert self.root.cell_count == len(cells) - outside_cells
        return outside_cells

    @classmethod
    def export_nodes_csv(cls, filename: str | Path, nodes: list[AtlasNode]) -> None:
        header = [
            "ID",
            "Abbreviation",
            "Name",
            "Total voxels",
            "Total volume (mm3)",
            "Total cells",
            "Cell density (mm3)",
        ]
        lines = [header]

        for node in nodes:
            line = [
                node.region_id,
                node.acronym,
                node.name,
                node.volume_voxels,
                node.volume_mm3,
                node.cell_count,
                node.cell_density_mm3,
            ]
            lines.append(list(map(str, line)))

        with open(filename, "w", newline="\n") as fh:
            writer = csv.writer(fh, delimiter=",")
            writer.writerows(lines)

    def export_to_vaa3d_format(self, filename: str | Path) -> None:
        root = self.root
        max_depth = max(d for _, d in root.pre_order_with_depth(depth=1))

        lines = []
        header = [
            "ID",
            "Abbreviation",
            "Leaf",
            "Parent ID",
            "Original ID",
            "Original abbreviation",
            "Total voxels",
            "Total volume (mm3)",
            "Total cells",
            "Cell density (mm3)",
        ]
        start_name_i = len(header)
        for i in range(max_depth):
            header.append(f"Region level {i + 1:02}")

        lines.append(header)
        for node, depth in root.pre_order_with_depth(depth=0):
            line = [
                "",
            ] * len(header)

            line[0] = str(node.region_id)
            line[1] = node.acronym
            line[2] = "FALSE" if node.children else "TRUE"
            line[3] = str(node.parent.region_id) if node.parent else "-1"

            line[4] = line[0]
            line[5] = line[1]
            if node.parent and node.parent.child_leaf is node:
                line[4] = str(node.parent.region_id)
                line[5] = node.parent.acronym

            line[6] = str(node.volume_voxels)
            line[7] = str(node.volume_mm3)
            line[8] = str(node.cell_count)
            line[9] = str(node.cell_density_mm3)

            line[start_name_i + depth] = node.name

            lines.append(line)

        with open(filename, "w", newline="\n") as fh:
            writer = csv.writer(fh, delimiter=",")
            writer.writerows(lines)

    def read_merge_tree(self, filename: str | Path) -> list[AtlasNode]:
        nodes = self.nodes_map

        with open(filename) as fh:
            reader = csv.reader(fh, delimiter=",")
            lines = [[c.strip() for c in row] for row in reader]

        header = lines[0]
        id_col = header.index("ID")
        acronym_col = header.index("Abbreviation")
        name_col = header.index("Description")
        children_col = max(id_col, max(acronym_col, name_col)) + 1

        merged_nodes = []
        for row in lines[1:]:
            node = AtlasNode(region_id=int(row[id_col]), acronym=row[acronym_col], name=row[name_col])
            items = [nodes[int(val)] for val in row[children_col:] if val]
            if not len(items):
                continue

            node.volume_voxels = sum(n.volume_voxels for n in items)
            node.volume_mm3 = sum(n.volume_mm3 for n in items)
            node.cell_count = sum(n.cell_count for n in items)

            merged_nodes.append(node)

        return merged_nodes
