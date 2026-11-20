"""Blender-independent tests for the custom molecule JSON data model."""
from __future__ import annotations

import ast
import collections
import copy
import math
import uuid
from pathlib import Path

ADDON = Path(__file__).resolve().parents[1] / "atomlab_addon.py"


def load_core():
    tree = ast.parse(ADDON.read_text(encoding="utf-8"), filename=str(ADDON))
    wanted_assignments = {
        "CUSTOM_MOLECULE_SCHEMA_VERSION", "_SUBSCRIPT_TRANSLATION",
        "_SUPERSCRIPT_TRANSLATION", "_LIBRARY_MARKUP_PATTERN",
    }
    wanted_functions = {
        "format_library_display_text", "custom_formula_from_symbols", "validate_custom_molecule_record",
        "validate_custom_library_payload", "merge_custom_library_records",
    }
    body = []
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id in wanted_assignments
            for target in node.targets
        ):
            body.append(node)
        elif isinstance(node, ast.FunctionDef) and node.name in wanted_functions:
            body.append(node)
    namespace = {
        "collections": collections, "math": math, "uuid": uuid, "re": __import__("re"),
        "ELEMENTS": [(1, "H", "Hydrogen", 1), (6, "C", "Carbon", 12),
                     (7, "N", "Nitrogen", 14), (8, "O", "Oxygen", 16)],
    }
    exec(compile(ast.Module(body=body, type_ignores=[]), str(ADDON), "exec"), namespace)
    return namespace


def water_record(name="Water"):
    return {
        "schema_version": 1,
        "id": uuid.uuid4().hex,
        "name": name,
        "formula": "H2O",
        "coordinate_unit": "blender_unit",
        "atoms": [
            {"element": "O", "charge": 0, "position": [0.0, 0.0, 0.0]},
            {"element": "H", "charge": 0, "position": [-0.7, 0.8, 0.0]},
            {"element": "H", "charge": 0, "position": [0.7, 0.8, 0.0]},
        ],
        "bonds": [{"a": 0, "b": 1, "order": 1}, {"a": 0, "b": 2, "order": 1}],
        "net_charge": 0,
        "note": "Test fixture",
    }


def run():
    c = load_core()
    assert c["custom_formula_from_symbols"](["O", "H", "H"]) == "H2O"
    assert c["custom_formula_from_symbols"](["O", "C", "H", "H", "H"]) == "CH3O"
    format_text = c["format_library_display_text"]
    assert format_text("H_2O") == "H₂O"
    assert format_text("SO_4^2-") == "SO₄²⁻"
    assert format_text("Fe^{3+}") == "Fe³⁺"
    assert format_text("follow_up note") == "follow_up note"

    record = water_record()
    assert c["validate_custom_molecule_record"](record) is record
    payload = {"schema_version": 1, "molecules": [record]}
    assert c["validate_custom_library_payload"](payload)["molecules"] == [record]
    assert c["merge_custom_library_records"]([], [record]) == [record]
    assert c["merge_custom_library_records"]([record], [record]) == [record]

    conflict = copy.deepcopy(record)
    conflict["id"] = uuid.uuid4().hex
    for label, mutation in [
        ("formula", lambda r: r.update(formula="H2O2")),
        ("net charge", lambda r: r.update(net_charge=1)),
        ("bond index", lambda r: r["bonds"][0].update(b=99)),
        ("duplicate bond", lambda r: r["bonds"].append(copy.deepcopy(r["bonds"][0]))),
        ("unsupported element", lambda r: r["atoms"][0].update(element="Xx")),
        ("nonfinite coordinate", lambda r: r["atoms"][0]["position"].__setitem__(0, float("nan"))),
        ("charge type", lambda r: r["atoms"][0].update(charge=True)),
    ]:
        bad = copy.deepcopy(record)
        mutation(bad)
        try:
            c["validate_custom_molecule_record"](bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"Validator accepted invalid {label}")

    try:
        c["merge_custom_library_records"]([record], [conflict])
    except ValueError:
        pass
    else:
        raise AssertionError("Name collision was overwritten without consent")
    replaced = c["merge_custom_library_records"]([record], [conflict], replace_conflicts=True)
    assert replaced == [conflict]

    try:
        c["validate_custom_library_payload"]({"schema_version": 99, "molecules": []})
    except ValueError:
        pass
    else:
        raise AssertionError("Unsupported library schema was accepted")

    print("Custom molecule schema tests passed: formula, charges, geometry, bonds, ids, merge/overwrite policy.")


if __name__ == "__main__":
    run()
