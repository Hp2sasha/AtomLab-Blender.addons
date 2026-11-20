# Test instructions

Run unit/regression checks from the repository root:

```sh
python tests/test_custom_molecule_schema.py
python tests/test_ionic_structures.py
```

Run Blender smoke tests with Blender 5.2 or later:

```sh
blender -b --factory-startup --python tests/blender_custom_molecule_smoke_test.py
blender -b --factory-startup --python tests/blender_ionic_smoke_test.py
```

The Blender custom-molecule test uses a temporary directory for its shared JSON library and disposable `.blend` project.
