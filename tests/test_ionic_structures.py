"""Blender-independent regression tests for AtomLab's ionic prototypes.

Run from the project directory with:
    python tests/test_ionic_structures.py

The addon imports bpy at module import time, so this harness extracts only its
pure-Python unit-cell data and lattice-expansion function from the AST.
"""
from __future__ import annotations

import ast
import collections
import math
from pathlib import Path

ADDON = Path(__file__).resolve().parents[1] / "atomlab_addon.py"

EXPECTED_COUNTS = {
    "NACL_ROCKSALT": {"Na": 32, "Cl": 32},
    "CSCL": {"Cs": 8, "Cl": 8},
    "ZNS_ZINCBLENDE": {"Zn": 32, "S": 32},
    "ZNS_WURTZITE": {"Zn": 16, "S": 16},
    "CAF2_FLUORITE": {"Ca": 32, "F": 64},
    "LI2O_ANTIFLUORITE": {"Li": 64, "O": 32},
    "TIO2_RUTILE": {"Ti": 16, "O": 32},
    "CATIO3_PEROVSKITE": {"Ca": 8, "Ti": 8, "O": 24},
    "NIAS_NICKEL_ARSENIDE": {"Ni": 16, "As": 16},
    "CDI2_LAYERED": {"Cd": 8, "I": 16},
    "MGAL2O4_SPINEL": {"Mg": 64, "Al": 128, "O": 256},
    "AL2O3_CORUNDUM": {"Al": 32, "O": 48},
}

# Expected nearest unlike-ion shell counts in the periodic bulk. Most models
# use a 3% tolerance; rutile's 4+2 and corundum's 3+3 distorted octahedra need
# data-specific tolerances. Perovskite and spinel oxygen have distinct shells.
EXPECTED_NEIGHBORS = {
    "NACL_ROCKSALT": {"Na": {"Cl": 6}, "Cl": {"Na": 6}},
    "CSCL": {"Cs": {"Cl": 8}, "Cl": {"Cs": 8}},
    "ZNS_ZINCBLENDE": {"Zn": {"S": 4}, "S": {"Zn": 4}},
    "ZNS_WURTZITE": {"Zn": {"S": 4}, "S": {"Zn": 4}},
    "CAF2_FLUORITE": {"Ca": {"F": 8}, "F": {"Ca": 4}},
    "LI2O_ANTIFLUORITE": {"Li": {"O": 4}, "O": {"Li": 8}},
    "TIO2_RUTILE": {"Ti": {"O": 6}, "O": {"Ti": 3}},
    "CATIO3_PEROVSKITE": {"Ca": {"O": 12}, "Ti": {"O": 6}, "O": {"Ti": 2, "Ca": 4}},
    "NIAS_NICKEL_ARSENIDE": {"Ni": {"As": 6}, "As": {"Ni": 6}},
    "CDI2_LAYERED": {"Cd": {"I": 6}, "I": {"Cd": 3}},
    "MGAL2O4_SPINEL": {"Mg": {"O": 4}, "Al": {"O": 6}, "O": {"Mg": 1, "Al": 3}},
    "AL2O3_CORUNDUM": {"Al": {"O": 6}, "O": {"Al": 4}},
}


def load_pure_addon_bits():
    tree = ast.parse(ADDON.read_text(encoding="utf-8"), filename=str(ADDON))
    wanted_assignments = {
        "_FCC_POSITIONS", "_TETRA_POSITIONS", "_SPINEL_MG_POSITIONS",
        "_SPINEL_AL_POSITIONS", "_SPINEL_O_POSITIONS", "IONIC_STRUCTURE_DATA",
        "IONIC_STRUCTURE_ITEMS",
    }
    wanted_functions = {"ionic_structure_sites"}
    body = []
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id in wanted_assignments
            for target in node.targets
        ):
            body.append(node)
        elif isinstance(node, ast.FunctionDef) and node.name in wanted_functions:
            body.append(node)
    namespace = {"math": math}
    exec(compile(ast.Module(body=body, type_ignores=[]), str(ADDON), "exec"), namespace)
    return namespace


def run():
    ns = load_pure_addon_bits()
    data = ns["IONIC_STRUCTURE_DATA"]
    items = ns["IONIC_STRUCTURE_ITEMS"]
    assert set(data) == set(EXPECTED_COUNTS), "prototype list changed; update regression expectations"
    assert {entry[0] for entry in items} == set(data), "UI selector/data mismatch"

    for key, record in data.items():
        sites = ns["ionic_structure_sites"](key, (2, 2, 2))
        counts = collections.Counter(element for element, _charge, _position in sites)
        assert dict(counts) == EXPECTED_COUNTS[key], (key, dict(counts))
        assert sum(charge for _element, charge, _position in sites) == 0, key
        coordinates = [tuple(round(v, 9) for v in pos) for _, _, pos in sites]
        assert len(coordinates) == len(set(coordinates)), f"duplicate sites: {key}"
        assert len(sites) == 8 * len(record["basis"]), key

        # Build a 3x3x3 supercell and inspect the central conventional cell,
        # avoiding artificial fragment-surface under-coordination.
        larger = ns["ionic_structure_sites"](key, (3, 3, 3))
        basis_size = len(record["basis"])
        center_cell = larger[13 * basis_size : 14 * basis_size]
        for element, charge, position in center_cell:
            for neighbor_element, expected_number in EXPECTED_NEIGHBORS[key][element].items():
                distances = sorted(
                    math.dist(position, other_position)
                    for other_element, other_charge, other_position in larger
                    if other_element == neighbor_element and charge * other_charge < 0
                )
                assert distances, (key, element, neighbor_element, "no unlike neighbors")
                threshold = distances[0] * (1.0 + record.get("coordination_tolerance", 0.03)) + 1e-9
                number = sum(distance <= threshold for distance in distances)
                assert number == expected_number, (
                    key, element, neighbor_element, number, expected_number, distances[:12]
                )

        print(f"PASS {key}: stoichiometry, neutral charge, unique sites, bulk coordination")

    print(f"All {len(EXPECTED_COUNTS)} ionic-prototype regression tests passed.")


if __name__ == "__main__":
    run()
