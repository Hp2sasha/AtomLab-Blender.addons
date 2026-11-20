"""In-Blender smoke test for AtomLab ionic structure creation.

Run only in a disposable/blank scene after enabling atomlab_addon.py. Each
preset replaces the previous AtomLab_Generated collection, so save your work
first. Set CONFIRM_DISPOSABLE_SCENE=True below, then run this script in
Blender's Scripting workspace. The final CaTiO3 fragment remains in the scene.
"""
import bpy
import collections

CONFIRM_DISPOSABLE_SCENE = False

EXPECTED = {
    "NACL_ROCKSALT": ("NaCl", {"Na": 32, "Cl": 32}),
    "CSCL": ("CsCl", {"Cs": 8, "Cl": 8}),
    "ZNS_ZINCBLENDE": ("ZnS", {"Zn": 32, "S": 32}),
    "ZNS_WURTZITE": ("ZnS", {"Zn": 16, "S": 16}),
    "CAF2_FLUORITE": ("CaF2", {"Ca": 32, "F": 64}),
    "LI2O_ANTIFLUORITE": ("Li2O", {"Li": 64, "O": 32}),
    "TIO2_RUTILE": ("TiO2", {"Ti": 16, "O": 32}),
    "CATIO3_PEROVSKITE": ("CaTiO3", {"Ca": 8, "Ti": 8, "O": 24}),
    "NIAS_NICKEL_ARSENIDE": ("NiAs", {"Ni": 16, "As": 16}),
    "CDI2_LAYERED": ("CdI2", {"Cd": 8, "I": 16}),
    "MGAL2O4_SPINEL": ("MgAl2O4", {"Mg": 64, "Al": 128, "O": 256}),
    "AL2O3_CORUNDUM": ("Al2O3", {"Al": 32, "O": 48}),
}


def run():
    if not CONFIRM_DISPOSABLE_SCENE:
        raise RuntimeError(
            "Safety stop: save your work, use a disposable/blank scene, then "
            "set CONFIRM_DISPOSABLE_SCENE=True in this script."
        )
    scene = bpy.context.scene
    if not hasattr(scene, "atomlab_ionic_structure"):
        raise RuntimeError("AtomLab is not enabled or its scene properties are missing.")

    lines = ["AtomLab ionic Blender smoke test", ""]
    for structure_id, (formula, expected_counts) in EXPECTED.items():
        scene.atomlab_ionic_structure = structure_id
        result = bpy.ops.atomlab.create_ionic_structure()
        if "FINISHED" not in result:
            raise AssertionError(f"Operator did not finish for {structure_id}: {result}")

        collection = bpy.data.collections.get("AtomLab_Generated")
        if collection is None:
            raise AssertionError(f"Generated collection missing for {structure_id}")
        objects = list(collection.objects)
        counts = collections.Counter(obj.get("element_symbol") for obj in objects)
        counts.pop(None, None)
        if dict(counts) != expected_counts:
            raise AssertionError(
                f"{structure_id}: expected {expected_counts}, got {dict(counts)}"
            )
        if any(obj.type != "MESH" for obj in objects):
            raise AssertionError(f"Unexpected non-mesh object (possible bond/rod) in {structure_id}")
        if any(obj.get("is_discrete_molecule") is not False for obj in objects):
            raise AssertionError(f"Molecular metadata incorrect for {structure_id}")
        if any("formal_charge" not in obj for obj in objects):
            raise AssertionError(f"Formal-charge metadata missing for {structure_id}")
        net_charge = sum(int(obj["formal_charge"]) for obj in objects)
        if net_charge != 0:
            raise AssertionError(f"{structure_id}: unexpected net charge {net_charge}")
        if scene.get("atomlab_ionic_structure_id") != structure_id:
            raise AssertionError(f"Scene structure metadata mismatch for {structure_id}")
        if scene.get("atomlab_formula") != formula:
            raise AssertionError(f"Scene formula metadata mismatch for {structure_id}")

        line = f"PASS {structure_id}: {len(objects)} ions, {dict(counts)}, net charge 0"
        print(line)
        lines.append(line)

    report = "\n".join(lines)
    text = bpy.data.texts.get("AtomLab_Ionic_Smoke_Test_Report")
    if text is None:
        text = bpy.data.texts.new("AtomLab_Ionic_Smoke_Test_Report")
    else:
        text.clear()
    text.write(report)
    print(f"All {len(EXPECTED)} Blender smoke tests passed. Report saved in Text datablock:", text.name)
    return report


run()
