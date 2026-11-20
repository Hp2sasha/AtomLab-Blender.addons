# AtomLab Blender Add-on

AtomLab is an educational Blender add-on for schematic atom models, ionic crystal fragments, and user-built molecules.

**Current version:** 0.16.0  
**Target:** Blender 5.2+

## Install

1. In Blender, open **Edit → Preferences → Add-ons**.
2. Use the top-right menu **Install from Disk…** and select `atomlab_addon.py`.
3. Enable **AtomLab** in the Add-ons list.
4. In the 3D View, press **N** and open the **AtomLab** sidebar tab.

## Custom molecule library

The add-on stores custom molecule records both in the current `.blend` file and in a shared local JSON library. JSON import/export is available in the panel. Records retain atom elements, coordinates, formal charges, bonds, bond orders, formula, net charge, notes, and a stable identifier. Name/ID conflicts are not replaced unless explicitly allowed.

In the library name and note fields, use lightweight plain-text notation:

- `H_2O` → H₂O
- `SO_4^2-` → SO₄²⁻
- `Fe^{3+}` → Fe³⁺

The preview and saved-molecule list render subscript/superscript characters; JSON retains the portable notation.

## Tests

- `tests/test_custom_molecule_schema.py`: Blender-independent validation, merge, and markup tests. Run with `python tests/test_custom_molecule_schema.py`.
- `tests/blender_custom_molecule_smoke_test.py`: Blender integration test for registration, storage, JSON import/export, and `.blend` reopen persistence. Run from the repository root with Blender 5.2+ in background mode, for example `blender -b --factory-startup --python tests/blender_custom_molecule_smoke_test.py`.
- `tests/test_ionic_structures.py`: regression checks for the 12 ionic structure prototypes.
- `tests/blender_ionic_smoke_test.py`: Blender integration test for ionic structures.

## Project status

Quantum-model review was approved; the persistent custom-molecule library is implemented and tested. The standalone application phase has not started. See `docs/` for the roadmap and project notes.

## License

No license has been selected yet. All rights reserved unless the repository owner adds a license.
