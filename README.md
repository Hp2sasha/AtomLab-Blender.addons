# AtomLab — Educational Blender Add-on

**Version:** 0.16.0  
**Target:** Blender 5.2 and later  
**Project type:** University Software Laboratory course project

> **Academic-project notice**  
> AtomLab is prepared solely as an educational project for presentation and assessment in a university Software Laboratory course. It is a course prototype—not a production application, certified scientific software, or an official product of any university or of the Blender Foundation. Its scientific visualizations are simplified and schematic; do not use them as measured, research-grade, or safety-critical data.

## Overview

AtomLab is a Blender add-on for creating and exploring educational 3D visualizations of atomic models, selected quantum-orbital probability shapes, extended ionic crystal fragments, and user-built molecules.

The add-on emphasizes clear educational limitations: idealized geometry, formal charges, illustrative scales, and schematic renderings are identified as such in the interface. Ionic structures are represented as extended crystal fragments; the add-on does not imply covalent rods between ions.

## Features

- Schematic atom and nucleus visualizations, including simplified quantum-orbital probability representations.
- A selection of extensible ionic-crystal prototypes represented as lattice fragments.
- Interactive molecule builder with atom elements, formal charges, schematic bonds, and bond orders.
- Persistent custom-molecule library stored both in the current `.blend` project and in a shared local JSON file.
- Versioned JSON import/export for moving custom-molecule records between installations.
- Validation of records and explicit protection against unintended name/ID replacement.
- Subscript/superscript notation in custom library names and notes:
  - `H_2O` displays as H₂O.
  - `SO_4^2-` displays as SO₄²⁻.
  - `Fe^{3+}` displays as Fe³⁺.

## Installation

1. In Blender, open **Edit → Preferences → Add-ons**.
2. Select **Install from Disk…** from the menu at the upper-right of the Add-ons section.
3. Choose `atomlab_addon.py` and enable **AtomLab** in the Add-ons list.
4. In the 3D Viewport, press **N** and select the **AtomLab** sidebar tab.

## Custom-molecule library

Use the **Interactive molecule builder** section to add atoms, assign formal charges, create bonds, and optionally optimize supported common geometries. In **My Molecule Library**, provide a name and optional note, then save the molecule.

Custom records preserve a stable identifier, name, derived formula, atom elements, 3D coordinates, per-atom formal charges, bond endpoints and orders, net charge, and note. Conflicting names or identifiers are rejected unless explicit replacement is enabled.

The project copy is embedded in the Blender Text datablock `AtomLab_Custom_Molecules.json`; save the `.blend` file to keep it. The shared library is stored in Blender's user configuration directory at `atomlab/custom_molecules.json`. Use the panel's JSON import/export controls to transfer records.

Subscript/superscript markup is stored as plain text for portability and rendered in the library preview and saved-molecule selector. Write `_` followed by the subscript, or `^` followed by the superscript; braces can delimit a group, for example `^{3+}`.

## Tests

Run the Blender-independent tests from this repository's root:

```sh
python tests/test_custom_molecule_schema.py
python tests/test_ionic_structures.py
```

Run the integration smoke tests with Blender 5.2 or later:

```sh
blender -b --factory-startup --python tests/blender_custom_molecule_smoke_test.py
blender -b --factory-startup --python tests/blender_ionic_smoke_test.py
```

The custom-molecule smoke test exercises registration, validation, JSON import/export, project embedding, and `.blend` reopen persistence. The ionic test script is intended for a disposable test scene; follow its in-file instructions before running it interactively.

## Repository layout

```text
atomlab_addon.py          Blender add-on source
README.md                 English project overview
README_FA.md              Persian project overview
docs/                     Roadmap, project plan, review checklist
tests/                    Unit, regression, and Blender smoke tests
```

## Project status

The quantum-model review was approved, and the persistent custom-molecule library has been implemented and smoke-tested on Blender 5.2.4. The standalone-application phase has not started. See `docs/AtomLab_Product_Roadmap.md` for the current roadmap and `docs/AtomLab_Project_Plan.md` for design context.

## License and reuse

No open-source license has been selected for this academic project. Unless a license is added by the project owner, publication on GitHub does not itself grant permission to reuse, redistribute, or modify the code. Blender is a separate project and is not affiliated with or endorsing AtomLab.
