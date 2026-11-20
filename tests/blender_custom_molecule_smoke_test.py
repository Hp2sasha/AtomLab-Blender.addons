"""Blender integration smoke test for custom-molecule persistence.

Run in a disposable scene after enabling the addon, or headlessly with the addon
module registered. It creates test objects and a temporary shared JSON library.
"""
import bpy
import json
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import atomlab_addon as addon


def run():
    if not hasattr(bpy.context.scene, "atomlab_custom_molecule_name"):
        addon.register()
    scene = bpy.context.scene
    temp_dir = tempfile.mkdtemp(prefix="atomlab_custom_test_")
    shared_path = os.path.join(temp_dir, "custom_molecules.json")
    blend_path = os.path.join(temp_dir, "project_copy.blend")
    addon.custom_library_user_path = lambda: shared_path

    # Fresh builder atoms for a bent water-like test graph.
    addon.clear_builder_collection(scene)
    atom_o = addon.add_builder_atom(bpy.context, 8, 0)
    atom_h1 = addon.add_builder_atom(bpy.context, 1, 0)
    atom_h2 = addon.add_builder_atom(bpy.context, 1, 0)
    atom_o.location = (0.0, 0.0, 0.0)
    atom_h1.location = (-0.7, 0.8, 0.0)
    atom_h2.location = (0.7, 0.8, 0.0)
    for atom_a, atom_b in ((atom_o, atom_h1), (atom_o, atom_h2)):
        bpy.ops.object.select_all(action='DESELECT')
        atom_a.select_set(True)
        atom_b.select_set(True)
        bpy.context.view_layer.objects.active = atom_a
        addon.create_builder_bond(bpy.context, 1)

    scene.atomlab_custom_molecule_name = "Smoke Test H_2O"
    scene.atomlab_custom_molecule_note = "Temporary Blender integration test"
    record = addon.custom_record_from_builder(
        bpy.context, scene.atomlab_custom_molecule_name, scene.atomlab_custom_molecule_note
    )
    assert record["formula"] == "H2O"
    assert len(record["atoms"]) == 3 and len(record["bonds"]) == 2
    assert record["net_charge"] == 0

    addon.save_custom_library_record_to_both(bpy.context, record)
    assert os.path.isfile(shared_path)
    with open(shared_path, encoding='utf-8') as handle:
        shared = addon.validate_custom_library_payload(json.load(handle))
    project = addon.read_project_custom_library()
    assert shared["molecules"] == [record]
    assert project["molecules"] == [record]

    addon.refresh_custom_molecule_catalog()
    enum_key = 'CUSTOM_' + record['id'].upper()
    assert enum_key in addon.CUSTOM_MOLECULES_BY_ID
    assert addon.format_library_display_text(record["name"]) == "Smoke Test H₂O"
    enum_label = next(item[1] for item in addon.CUSTOM_MOLECULE_ITEMS if item[0] == enum_key)
    assert enum_label == "Smoke Test H₂O"
    created = addon.create_saved_custom_molecule(bpy.context, record)
    assert created["id"] == record["id"]
    collection = next(c for c in bpy.data.collections if c.get("atomlab_custom_molecule_id") == record["id"])
    assert sum(1 for obj in collection.objects if obj.get("element_symbol")) == 3
    assert sum(1 for obj in collection.objects if obj.get("atomlab_custom_molecule_id") == record["id"] and "bond_order" in obj) == 2

    export_path = os.path.join(temp_dir, "roundtrip_export.json")
    export_result = bpy.ops.atomlab.export_custom_molecule_library(
        'EXEC_DEFAULT', filepath=export_path
    )
    assert 'FINISHED' in export_result
    with open(export_path, encoding='utf-8') as handle:
        exported = addon.validate_custom_library_payload(json.load(handle))
    assert len(exported["molecules"]) == 1
    extra = dict(record)
    extra["id"] = "4dcb683d9c6c42bfa032000000000001"
    extra["name"] = "Imported Test Water"
    import_path = os.path.join(temp_dir, "roundtrip_import.json")
    with open(import_path, 'w', encoding='utf-8') as handle:
        json.dump({"schema_version": 1, "molecules": [extra]}, handle)
    import_result = bpy.ops.atomlab.import_custom_molecule_library(
        'EXEC_DEFAULT', filepath=import_path
    )
    assert 'FINISHED' in import_result
    assert len(addon.read_project_custom_library()["molecules"]) == 2
    assert len(addon.read_global_custom_library()["molecules"]) == 2

    bpy.ops.wm.save_as_mainfile(filepath=blend_path)
    # Reopen the saved project and verify the embedded JSON plus rebuilt scene data survive.
    bpy.ops.wm.open_mainfile(filepath=blend_path)
    restored = addon.read_project_custom_library()["molecules"]
    assert len(restored) == 2 and restored[0]["id"] == record["id"]
    restored_record = next(item for item in restored if item["id"] == record["id"])
    assert restored_record["atoms"] == record["atoms"]
    restored_collection = next(
        c for c in bpy.data.collections if c.get("atomlab_custom_molecule_id") == record["id"]
    )
    restored_atoms = [
        obj for obj in restored_collection.objects
        if obj.get("atomlab_custom_molecule_id") == record["id"] and obj.get("element_symbol")
    ]
    restored_bonds = [
        obj for obj in restored_collection.objects
        if obj.get("atomlab_custom_molecule_id") == record["id"] and "bond_order" in obj
    ]
    assert len(restored_atoms) == 3 and len(restored_bonds) == 2
    for obj in restored_atoms:
        idx = int(obj["atomlab_custom_atom_index"])
        assert obj["element_symbol"] == record["atoms"][idx]["element"]
        assert int(obj["formal_charge"]) == record["atoms"][idx]["charge"]
        assert all(abs(float(a) - float(b)) < 1e-5 for a, b in zip(obj.location, record["atoms"][idx]["position"]))
    assert {obj["bond_order"] for obj in restored_bonds} == {1}
    assert len(addon.read_global_custom_library()["molecules"]) == 2
    addon.unregister()
    addon.register()
    addon.unregister()
    print("CUSTOM_MOLECULE_SMOKE_PASS", shared_path, blend_path)
    return temp_dir, blend_path


if __name__ == "__main__":
    run()
