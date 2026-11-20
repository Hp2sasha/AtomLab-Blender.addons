bl_info = {
    "name": "AtomLab — Educational Atom Models",
    "author": "AtomLab Project",
    "version": (0, 16, 0),
    "blender": (5, 2, 0),
    "location": "View3D > Sidebar > AtomLab",
    "description": "Educational, explicitly schematic atom visualizer",
    "category": "3D View",
}

import bpy
import collections
import json
import math
import os
import random
import re
import tempfile
import uuid
from mathutils import Vector
from bpy.props import BoolProperty, EnumProperty, IntProperty, StringProperty


# Neutral atoms Z=1..118. The listed mass is a representative isotope/mass number;
# for radioactive elements it is a commonly tabulated nuclide, not an abundance claim.
ELEMENT_ROWS = """
1 H Hydrogen 1
2 He Helium 4
3 Li Lithium 7
4 Be Beryllium 9
5 B Boron 11
6 C Carbon 12
7 N Nitrogen 14
8 O Oxygen 16
9 F Fluorine 19
10 Ne Neon 20
11 Na Sodium 23
12 Mg Magnesium 24
13 Al Aluminium 27
14 Si Silicon 28
15 P Phosphorus 31
16 S Sulfur 32
17 Cl Chlorine 35
18 Ar Argon 40
19 K Potassium 39
20 Ca Calcium 40
21 Sc Scandium 45
22 Ti Titanium 48
23 V Vanadium 51
24 Cr Chromium 52
25 Mn Manganese 55
26 Fe Iron 56
27 Co Cobalt 59
28 Ni Nickel 58
29 Cu Copper 63
30 Zn Zinc 64
31 Ga Gallium 69
32 Ge Germanium 74
33 As Arsenic 75
34 Se Selenium 80
35 Br Bromine 79
36 Kr Krypton 84
37 Rb Rubidium 85
38 Sr Strontium 88
39 Y Yttrium 89
40 Zr Zirconium 90
41 Nb Niobium 93
42 Mo Molybdenum 98
43 Tc Technetium 98
44 Ru Ruthenium 102
45 Rh Rhodium 103
46 Pd Palladium 106
47 Ag Silver 107
48 Cd Cadmium 114
49 In Indium 115
50 Sn Tin 120
51 Sb Antimony 121
52 Te Tellurium 130
53 I Iodine 127
54 Xe Xenon 132
55 Cs Caesium 133
56 Ba Barium 138
57 La Lanthanum 139
58 Ce Cerium 140
59 Pr Praseodymium 141
60 Nd Neodymium 142
61 Pm Promethium 145
62 Sm Samarium 152
63 Eu Europium 153
64 Gd Gadolinium 158
65 Tb Terbium 159
66 Dy Dysprosium 164
67 Ho Holmium 165
68 Er Erbium 166
69 Tm Thulium 169
70 Yb Ytterbium 174
71 Lu Lutetium 175
72 Hf Hafnium 180
73 Ta Tantalum 181
74 W Tungsten 184
75 Re Rhenium 187
76 Os Osmium 192
77 Ir Iridium 193
78 Pt Platinum 195
79 Au Gold 197
80 Hg Mercury 202
81 Tl Thallium 205
82 Pb Lead 208
83 Bi Bismuth 209
84 Po Polonium 209
85 At Astatine 210
86 Rn Radon 222
87 Fr Francium 223
88 Ra Radium 226
89 Ac Actinium 227
90 Th Thorium 232
91 Pa Protactinium 231
92 U Uranium 238
93 Np Neptunium 237
94 Pu Plutonium 244
95 Am Americium 243
96 Cm Curium 247
97 Bk Berkelium 247
98 Cf Californium 251
99 Es Einsteinium 252
100 Fm Fermium 257
101 Md Mendelevium 258
102 No Nobelium 259
103 Lr Lawrencium 266
104 Rf Rutherfordium 267
105 Db Dubnium 268
106 Sg Seaborgium 269
107 Bh Bohrium 270
108 Hs Hassium 277
109 Mt Meitnerium 278
110 Ds Darmstadtium 281
111 Rg Roentgenium 282
112 Cn Copernicium 285
113 Nh Nihonium 286
114 Fl Flerovium 289
115 Mc Moscovium 290
116 Lv Livermorium 293
117 Ts Tennessine 294
118 Og Oganesson 294
"""
ELEMENTS = [
    (int(row.split()[0]), row.split()[1], row.split()[2], int(row.split()[3]))
    for row in ELEMENT_ROWS.strip().splitlines()
]

# Madelung (Aufbau) filling order through 7p, plus 8s for simple anion extension;
# selected neutral-atom exceptions are corrected below. Not a precision calculation.
AUFBAU_ORBITALS = (
    (1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2), (3, 1, 6),
    (4, 0, 2), (3, 2, 10), (4, 1, 6), (5, 0, 2), (4, 2, 10),
    (5, 1, 6), (6, 0, 2), (4, 3, 14), (5, 2, 10), (6, 1, 6),
    (7, 0, 2), (5, 3, 14), (6, 2, 10), (7, 1, 6), (8, 0, 2),
)
# (Z, subshell changes) to match compiled neutral-atom ground states.
# Sources: NIST electronic-configuration compilation (H-U) and NIST ASD ground-shell
# data for the later actinides. See AtomLab_Project_Plan.md for source URLs.
CONFIGURATION_ADJUSTMENTS = {
    24: {(4, 0): -1, (3, 2): 1}, 29: {(4, 0): -1, (3, 2): 1},
    41: {(5, 0): -1, (4, 2): 1}, 42: {(5, 0): -1, (4, 2): 1},
    44: {(5, 0): -1, (4, 2): 1}, 45: {(5, 0): -1, (4, 2): 1},
    46: {(5, 0): -2, (4, 2): 2}, 47: {(5, 0): -1, (4, 2): 1},
    57: {(4, 3): -1, (5, 2): 1}, 58: {(4, 3): -1, (5, 2): 1},
    64: {(4, 3): -1, (5, 2): 1}, 78: {(6, 0): -1, (5, 2): 1},
    79: {(6, 0): -1, (5, 2): 1}, 89: {(5, 3): -1, (6, 2): 1},
    90: {(5, 3): -2, (6, 2): 2}, 91: {(5, 3): -1, (6, 2): 1},
    92: {(5, 3): -1, (6, 2): 1}, 93: {(5, 3): -1, (6, 2): 1},
    96: {(5, 3): -1, (6, 2): 1},
    # NIST ASD identifies the Lr neutral ground shell as 5f14 7s2 7p1.
    103: {(6, 2): -1, (7, 1): 1},
}

ORBITAL_SHAPES = {
    0: [('S', None, 's', 'Spherical orbital')],
    1: [('PX', 'X', 'p_x', 'p orbital oriented along X'),
        ('PY', 'Y', 'p_y', 'p orbital oriented along Y'),
        ('PZ', 'Z', 'p_z', 'p orbital oriented along Z')],
    2: [('DXY', 'DXY', 'd_xy', 'real d_xy orbital shape'),
        ('DXZ', 'DXZ', 'd_xz', 'real d_xz orbital shape'),
        ('DYZ', 'DYZ', 'd_yz', 'real d_yz orbital shape'),
        ('DX2Y2', 'DX2Y2', 'd_x2-y2', 'real d_x2-y2 orbital shape'),
        ('DZ2', 'DZ2', 'd_z2', 'real d_z2 orbital shape')],
    3: [('FZ3', 'FZ3', 'f_z(5z^2-3r^2)', 'real f_z(5z^2-3r^2) orbital shape'),
        ('FXZ2', 'FXZ2', 'f_x(5z^2-r^2)', 'real f_x(5z^2-r^2) orbital shape'),
        ('FYZ2', 'FYZ2', 'f_y(5z^2-r^2)', 'real f_y(5z^2-r^2) orbital shape'),
        ('FZXY', 'FZXY', 'f_z(x^2-y^2)', 'real f_z(x^2-y^2) orbital shape'),
        ('FXYZ', 'FXYZ', 'f_xyz', 'real f_xyz orbital shape'),
        ('FXX23Y2', 'FXX23Y2', 'f_x(x^2-3y^2)', 'real f_x(x^2-3y^2) orbital shape'),
        ('FY3X2Y2', 'FY3X2Y2', 'f_y(3x^2-y^2)', 'real f_y(3x^2-y^2) orbital shape')],
}
ORBITAL_ITEMS = []
ORBITAL_DEFINITIONS = {}
for _n, _l, _capacity in AUFBAU_ORBITALS:
    for _suffix, _orientation, _shape, _description in ORBITAL_SHAPES[_l]:
        _orbital_id = f"{_n}{_suffix}"
        _orbital_name = f"{_n}{_shape}"
        ORBITAL_ITEMS.append((
            _orbital_id,
            _orbital_name,
            f"Representative {_n}{_description}",
        ))
        ORBITAL_DEFINITIONS[_orbital_id] = (
            _n, _l, _orientation, _orbital_name
        )


def _water_coordinates():
    half = math.radians(104.5 / 2.0)
    bond = 0.9584
    return [
        ('O', (0.0, 0.0, 0.0)),
        ('H', (-bond * math.sin(half), bond * math.cos(half), 0.0)),
        ('H', ( bond * math.sin(half), bond * math.cos(half), 0.0)),
    ]


def _methane_coordinates():
    bond = 1.09
    corners = ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))
    return [('C', (0.0, 0.0, 0.0))] + [
        ('H', tuple(bond * c / math.sqrt(3.0) for c in corner))
        for corner in corners
    ]


def _ammonia_coordinates():
    bond = 1.012
    angle = math.radians(107.8)
    cos_polar = math.sqrt((1.0 + 2.0 * math.cos(angle)) / 3.0)
    sin_polar = math.sqrt(1.0 - cos_polar * cos_polar)
    atoms = [('N', (0.0, 0.0, 0.0))]
    for azimuth in (0.0, 2.0 * math.pi / 3.0, 4.0 * math.pi / 3.0):
        atoms.append(('H', (
            bond * sin_polar * math.cos(azimuth),
            bond * sin_polar * math.sin(azimuth),
            bond * cos_polar,
        )))
    return atoms


MOLECULE_DATA = {
    'H2': {
        'label': 'Hydrogen (H₂)', 'formula': 'H2', 'geometry': 'linear',
        'atoms': [('H', (-0.3707, 0.0, 0.0)), ('H', (0.3707, 0.0, 0.0))],
        'bonds': [(0, 1, 1)], 'bond_length_A': 0.7414,
    },
    'O2': {
        'label': 'Oxygen (O₂)', 'formula': 'O2', 'geometry': 'linear',
        'atoms': [('O', (-0.6038, 0.0, 0.0)), ('O', (0.6038, 0.0, 0.0))],
        'bonds': [(0, 1, 2)], 'bond_length_A': 1.2075,
    },
    'H2O': {
        'label': 'Water (H₂O)', 'formula': 'H2O', 'geometry': 'bent, 104.5°',
        'atoms': _water_coordinates(), 'bonds': [(0, 1, 1), (0, 2, 1)],
        'bond_length_A': 0.9584,
    },
    'CO2': {
        'label': 'Carbon dioxide (CO₂)', 'formula': 'CO2', 'geometry': 'linear',
        'atoms': [('O', (-1.16, 0.0, 0.0)), ('C', (0.0, 0.0, 0.0)),
                  ('O', (1.16, 0.0, 0.0))],
        'bonds': [(0, 1, 2), (1, 2, 2)], 'bond_length_A': 1.16,
    },
    'CH4': {
        'label': 'Methane (CH₄)', 'formula': 'CH4', 'geometry': 'tetrahedral, 109.47°',
        'atoms': _methane_coordinates(),
        'bonds': [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)],
        'bond_length_A': 1.09,
    },
    'NH3': {
        'label': 'Ammonia (NH₃)', 'formula': 'NH3', 'geometry': 'trigonal pyramidal, 107.8°',
        'atoms': _ammonia_coordinates(),
        'bonds': [(0, 1, 1), (0, 2, 1), (0, 3, 1)],
        'bond_length_A': 1.012,
    },
}
MOLECULE_ITEMS = [
    (key, data['label'], f"{data['formula']} — {data['geometry']}")
    for key, data in MOLECULE_DATA.items()
]
CUSTOM_MOLECULE_SCHEMA_VERSION = 1
CUSTOM_MOLECULE_TEXT_NAME = "AtomLab_Custom_Molecules.json"
CUSTOM_MOLECULE_NONE = "__ATOMLAB_NO_CUSTOM_MOLECULE__"
CUSTOM_MOLECULE_ITEMS = [(CUSTOM_MOLECULE_NONE, "No saved molecules", "Save a molecule from the Builder first")]
CUSTOM_MOLECULES_BY_ID = {}
_SUBSCRIPT_TRANSLATION = str.maketrans("0123456789+-=()", "₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎")
_SUPERSCRIPT_TRANSLATION = str.maketrans("0123456789+-=()", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾")
_LIBRARY_MARKUP_PATTERN = re.compile(r'([_^])(?:\{([^{}]+)\}|([0-9+\-=()]+))')


def format_library_display_text(value):
    """Render lightweight _subscript and ^superscript notation for Blender labels.

    Examples: H_2O -> H₂O, SO_4^2- -> SO₄²⁻, Fe^{3+} -> Fe³⁺.
    The stored value remains plain text with markup for portable JSON round-trips.
    """
    translations = {'_': _SUBSCRIPT_TRANSLATION, '^': _SUPERSCRIPT_TRANSLATION}

    def replace(match):
        marker = match.group(1)
        content = match.group(2) if match.group(2) is not None else match.group(3)
        return content.translate(translations[marker])

    return _LIBRARY_MARKUP_PATTERN.sub(replace, str(value))


def custom_formula_from_symbols(symbols):
    """Return an uncomplicated Hill-style molecular formula string."""
    counts = collections.Counter(symbols)
    order = []
    if 'C' in counts:
        order.append('C')
        if 'H' in counts:
            order.append('H')
    order.extend(sorted(symbol for symbol in counts if symbol not in order))
    return ''.join(symbol + (str(counts[symbol]) if counts[symbol] != 1 else '')
                   for symbol in order)


def validate_custom_molecule_record(record):
    """Validate a versioned, bpy-independent custom molecule record."""
    if not isinstance(record, dict):
        raise ValueError("Molecule entry must be a JSON object")
    if record.get('schema_version') != CUSTOM_MOLECULE_SCHEMA_VERSION:
        raise ValueError("Unsupported molecule schema_version")
    if not isinstance(record.get('id'), str) or not record['id'].strip() or len(record['id']) > 128:
        raise ValueError("Molecule id must be a non-empty string of at most 128 characters")
    try:
        uuid.UUID(record['id'])
    except (ValueError, AttributeError):
        raise ValueError("Molecule id must be a UUID (generated automatically by AtomLab).")
    if not isinstance(record.get('name'), str) or not record['name'].strip():
        raise ValueError("Molecule name is required")
    atoms = record.get('atoms')
    bonds = record.get('bonds')
    if not isinstance(atoms, list) or not atoms or len(atoms) > 512:
        raise ValueError("Molecule must contain between 1 and 512 atoms")
    if not isinstance(bonds, list):
        raise ValueError("Molecule bonds must be a list")
    symbols = []
    total_charge = 0
    for index, atom in enumerate(atoms):
        if not isinstance(atom, dict):
            raise ValueError(f"Atom {index + 1} must be a JSON object")
        symbol = atom.get('element')
        if symbol not in {item[1] for item in ELEMENTS}:
            raise ValueError(f"Atom {index + 1} has unsupported element {symbol!r}")
        charge = atom.get('charge')
        if type(charge) is not int or not -8 <= charge <= 8:
            raise ValueError(f"Atom {index + 1} formal charge must be an integer from -8 to +8")
        position = atom.get('position')
        if not isinstance(position, list) or len(position) != 3:
            raise ValueError(f"Atom {index + 1} position must contain three coordinates")
        if any(type(value) not in (int, float) or not math.isfinite(value) for value in position):
            raise ValueError(f"Atom {index + 1} position contains a non-finite coordinate")
        symbols.append(symbol)
        total_charge += charge
    expected_formula = custom_formula_from_symbols(symbols)
    if record.get('formula') != expected_formula:
        raise ValueError(f"Formula mismatch; expected {expected_formula}")
    if type(record.get('net_charge')) is not int or record['net_charge'] != total_charge:
        raise ValueError(f"Net charge mismatch; expected {total_charge:+d}")
    if record.get('coordinate_unit') != 'blender_unit':
        raise ValueError("coordinate_unit must be 'blender_unit'")
    seen_pairs = set()
    for index, bond in enumerate(bonds):
        if not isinstance(bond, dict):
            raise ValueError(f"Bond {index + 1} must be a JSON object")
        first, second, order = bond.get('a'), bond.get('b'), bond.get('order')
        if type(first) is not int or type(second) is not int or first == second:
            raise ValueError(f"Bond {index + 1} must reference two distinct atom indices")
        if not (0 <= first < len(atoms) and 0 <= second < len(atoms)):
            raise ValueError(f"Bond {index + 1} references an atom outside the molecule")
        if type(order) is not int or order not in (1, 2, 3):
            raise ValueError(f"Bond {index + 1} order must be 1, 2, or 3")
        pair = tuple(sorted((first, second)))
        if pair in seen_pairs:
            raise ValueError(f"Duplicate bond between atoms {first} and {second}")
        seen_pairs.add(pair)
    if not isinstance(record.get('note', ''), str):
        raise ValueError("Molecule note must be text")
    return record


def validate_custom_library_payload(payload):
    if not isinstance(payload, dict) or payload.get('schema_version') != CUSTOM_MOLECULE_SCHEMA_VERSION:
        raise ValueError("Unsupported custom-library file or schema_version")
    molecules = payload.get('molecules')
    if not isinstance(molecules, list):
        raise ValueError("Library 'molecules' must be a list")
    seen_ids = set()
    validated = []
    for record in molecules:
        validate_custom_molecule_record(record)
        if record['id'] in seen_ids:
            raise ValueError(f"Duplicate molecule id in library: {record['id']}")
        seen_ids.add(record['id'])
        validated.append(record)
    return {'schema_version': CUSTOM_MOLECULE_SCHEMA_VERSION, 'molecules': validated}


def merge_custom_library_records(existing, incoming, replace_conflicts=False):
    """Merge records by id/name; never overwrite a conflicting record implicitly."""
    result = list(existing)
    for candidate in incoming:
        validate_custom_molecule_record(candidate)
        by_id = [i for i, old in enumerate(result) if old['id'] == candidate['id']]
        by_name = [i for i, old in enumerate(result)
                   if old['name'].strip().casefold() == candidate['name'].strip().casefold()]
        conflicts = sorted(set(by_id + by_name))
        if conflicts and not replace_conflicts:
            identical = len(conflicts) == 1 and result[conflicts[0]] == candidate
            if identical:
                continue
            raise ValueError(
                f"A molecule named {candidate['name']!r} or with the same id already exists; "
                "enable explicit replacement to overwrite it."
            )
        for index in reversed(conflicts):
            result.pop(index)
        result.append(candidate)
    return result


_FCC_POSITIONS = (
    (0.0, 0.0, 0.0), (0.0, 0.5, 0.5),
    (0.5, 0.0, 0.5), (0.5, 0.5, 0.0),
)
_TETRA_POSITIONS = tuple(
    (x, y, z) for x in (0.25, 0.75) for y in (0.25, 0.75) for z in (0.25, 0.75)
)
# MgAl2O4 spinel conventional-cell Wyckoff orbits, derived from the ordered
# Fd-3m prototype (Mg 8a, Al 16d, O 32e; oxygen x=0.26171).
_SPINEL_MG_POSITIONS = (
    (0.125, 0.125, 0.125), (0.125, 0.625, 0.625),
    (0.375, 0.375, 0.875), (0.375, 0.875, 0.375),
    (0.625, 0.125, 0.625), (0.625, 0.625, 0.125),
    (0.875, 0.375, 0.375), (0.875, 0.875, 0.875),
)
_SPINEL_AL_POSITIONS = (
    (0.0, 0.0, 0.5), (0.0, 0.25, 0.75), (0.0, 0.5, 0.0), (0.0, 0.75, 0.25),
    (0.25, 0.0, 0.75), (0.25, 0.25, 0.5), (0.25, 0.5, 0.25), (0.25, 0.75, 0.0),
    (0.5, 0.0, 0.0), (0.5, 0.25, 0.25), (0.5, 0.5, 0.5), (0.5, 0.75, 0.75),
    (0.75, 0.0, 0.25), (0.75, 0.25, 0.0), (0.75, 0.5, 0.75), (0.75, 0.75, 0.5),
)
_SPINEL_O_POSITIONS = (
    (0.01171, 0.01171, 0.73829), (0.01171, 0.23829, 0.51171),
    (0.01171, 0.51171, 0.23829), (0.01171, 0.73829, 0.01171),
    (0.23829, 0.01171, 0.51171), (0.23829, 0.23829, 0.73829),
    (0.23829, 0.51171, 0.01171), (0.23829, 0.73829, 0.23829),
    (0.26171, 0.26171, 0.26171), (0.26171, 0.48829, 0.48829),
    (0.26171, 0.76171, 0.76171), (0.26171, 0.98829, 0.98829),
    (0.48829, 0.26171, 0.48829), (0.48829, 0.48829, 0.26171),
    (0.48829, 0.76171, 0.98829), (0.48829, 0.98829, 0.76171),
    (0.51171, 0.01171, 0.23829), (0.51171, 0.23829, 0.01171),
    (0.51171, 0.51171, 0.73829), (0.51171, 0.73829, 0.51171),
    (0.73829, 0.01171, 0.01171), (0.73829, 0.23829, 0.23829),
    (0.73829, 0.51171, 0.51171), (0.73829, 0.73829, 0.73829),
    (0.76171, 0.26171, 0.76171), (0.76171, 0.48829, 0.98829),
    (0.76171, 0.76171, 0.26171), (0.76171, 0.98829, 0.48829),
    (0.98829, 0.26171, 0.98829), (0.98829, 0.48829, 0.76171),
    (0.98829, 0.76171, 0.48829), (0.98829, 0.98829, 0.26171),
)
IONIC_STRUCTURE_DATA = {
    'NACL_ROCKSALT': {
        'label': 'NaCl crystal fragment (rock salt)', 'formula': 'NaCl',
        'prototype': 'rock salt, cubic Fm-3m', 'formula_units_per_cell': 4,
        'coordination': 'Na: 6; Cl: 6',
        'cell': ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        'basis': [('Cl', -1, p) for p in _FCC_POSITIONS] + [
            ('Na', 1, p) for p in ((0.5, 0, 0), (0, 0.5, 0), (0, 0, 0.5), (0.5, 0.5, 0.5))
        ],
        'note': 'Extended ionic lattice; no discrete NaCl molecules or covalent rods.',
    },
    'CSCL': {
        'label': 'CsCl crystal fragment (CsCl type)', 'formula': 'CsCl',
        'prototype': 'CsCl type, cubic Pm-3m', 'formula_units_per_cell': 1,
        'coordination': 'Cs: 8; Cl: 8',
        'cell': ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        'basis': [('Cs', 1, (0, 0, 0)), ('Cl', -1, (0.5, 0.5, 0.5))],
        'note': 'Alternating 1:1 ionic crystal; finite surfaces have reduced coordination.',
    },
    'ZNS_ZINCBLENDE': {
        'label': 'ZnS crystal fragment (zinc blende)', 'formula': 'ZnS',
        'prototype': 'zinc blende, cubic F-43m', 'formula_units_per_cell': 4,
        'coordination': 'Zn: 4; S: 4',
        'cell': ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        'basis': [('S', -2, p) for p in _FCC_POSITIONS] + [
            ('Zn', 2, tuple((p[i] + 0.25) % 1.0 for i in range(3)))
            for p in _FCC_POSITIONS
        ],
        'note': 'Tetrahedral ZnS has substantial covalent character; ±2 values are formal oxidation-state labels.'
    },
    'ZNS_WURTZITE': {
        'label': 'ZnS crystal fragment (wurtzite)', 'formula': 'ZnS',
        'prototype': 'wurtzite, hexagonal P63mc', 'formula_units_per_cell': 2,
        'coordination': 'Zn: 4; S: 4',
        'cell': ((1, 0, 0), (-0.5, math.sqrt(3) / 2, 0), (0, 0, math.sqrt(8 / 3))),
        'basis': [
            ('S', -2, (1/3, 2/3, 0)), ('S', -2, (2/3, 1/3, 0.5)),
            ('Zn', 2, (1/3, 2/3, 0.375)), ('Zn', 2, (2/3, 1/3, 0.875)),
        ],
        'note': 'Idealized tetrahedral wurtzite ZnS; mixed bonding, normalized cell geometry, and formal ±2 labels.'
    },
    'CAF2_FLUORITE': {
        'label': 'CaF2 crystal fragment (fluorite)', 'formula': 'CaF2',
        'prototype': 'fluorite, cubic Fm-3m', 'formula_units_per_cell': 4,
        'coordination': 'Ca: 8; F: 4',
        'cell': ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        'basis': [('Ca', 2, p) for p in _FCC_POSITIONS] + [
            ('F', -1, p) for p in _TETRA_POSITIONS
        ],
        'note': 'CaF2 fluorite topology; formal charges +2 and -1 give charge neutrality.',
    },
    'LI2O_ANTIFLUORITE': {
        'label': 'Li2O crystal fragment (antifluorite)', 'formula': 'Li2O',
        'prototype': 'antifluorite, cubic Fm-3m', 'formula_units_per_cell': 4,
        'coordination': 'Li: 4; O: 8',
        'cell': ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        'basis': [('O', -2, p) for p in _FCC_POSITIONS] + [
            ('Li', 1, p) for p in _TETRA_POSITIONS
        ],
        'note': 'Li2O antifluorite topology; formal charges +1 and -2 give charge neutrality.',
    },
    'AL2O3_CORUNDUM': {
        'label': 'Al2O3 crystal fragment (corundum)', 'formula': 'Al2O3',
        'prototype': 'corundum, trigonal R-3c (primitive rhombohedral cell)',
        'formula_units_per_cell': 2,
        'coordination': 'Al: 6 octahedral; O: 4 Al',
        'coordination_tolerance': 0.07,
        'cell': (
            (1, 0, 0),
            (0.5692795, 0.8221440, 0),
            (0.5692795, 0.2982450, 0.7661402),
        ),
        'basis': [
            ('Al', 3, (0.3520, 0.3520, 0.3520)),
            ('Al', 3, (0.8520, 0.8520, 0.8520)),
            ('Al', 3, (0.6480, 0.6480, 0.6480)),
            ('Al', 3, (0.1480, 0.1480, 0.1480)),
            ('O', -2, (0.5563, 0.9437, 0.2500)),
            ('O', -2, (0.2500, 0.5563, 0.9437)),
            ('O', -2, (0.9437, 0.2500, 0.5563)),
            ('O', -2, (0.4437, 0.0563, 0.7500)),
            ('O', -2, (0.7500, 0.4437, 0.0563)),
            ('O', -2, (0.0563, 0.7500, 0.4437)),
        ],
        'note': 'Idealized alpha-Al2O3 corundum using the primitive rhombohedral cell and reported internal coordinates; formal charges are educational.'
    },
    'NIAS_NICKEL_ARSENIDE': {
        'label': 'NiAs crystal fragment (nickel arsenide)', 'formula': 'NiAs',
        'prototype': 'NiAs type, hexagonal P63/mmc', 'formula_units_per_cell': 2,
        'coordination': 'Ni: 6 octahedral; As: 6 trigonal-prismatic',
        'cell': ((1, 0, 0), (-0.5, math.sqrt(3) / 2, 0), (0, 0, math.sqrt(8 / 3))),
        'basis': [
            ('Ni', 3, (0, 0, 0)), ('Ni', 3, (0, 0, 0.5)),
            ('As', -3, (1/3, 2/3, 0.25)), ('As', -3, (2/3, 1/3, 0.75)),
        ],
        'note': 'Hexagonal NiAs prototype. NiAs has substantial metallic/covalent character; ±3 are formal teaching labels, not literal point charges.'
    },
    'CDI2_LAYERED': {
        'label': 'CdI2 crystal fragment (layered CdI2 type)', 'formula': 'CdI2',
        'prototype': 'CdI2 layered, trigonal P-3m1', 'formula_units_per_cell': 1,
        'coordination': 'Cd: 6 octahedral; I: 3 (layered)',
        'cell': ((1, 0, 0), (-0.5, math.sqrt(3) / 2, 0), (0, 0, 1.8)),
        'basis': [
            ('Cd', 2, (0, 0, 0)),
            ('I', -1, (1/3, 2/3, 0.25)), ('I', -1, (2/3, 1/3, 0.75)),
        ],
        'note': 'Idealized layered CdI2 topology with edge-sharing CdI6 octahedra; formal charges +2/-1. Interlayer cohesion is not drawn.'
    },
    'MGAL2O4_SPINEL': {
        'label': 'MgAl2O4 crystal fragment (spinel)', 'formula': 'MgAl2O4',
        'prototype': 'normal spinel, cubic Fd-3m', 'formula_units_per_cell': 8,
        'coordination': 'Mg: 4 tetrahedral; Al: 6 octahedral; O: 1 Mg + 3 Al',
        'cell': ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        'basis': (
            [('Mg', 2, p) for p in _SPINEL_MG_POSITIONS]
            + [('Al', 3, p) for p in _SPINEL_AL_POSITIONS]
            + [('O', -2, p) for p in _SPINEL_O_POSITIONS]
        ),
        'note': 'Ordered normal-spinel topology; oxygen internal coordinate x=0.26171. Formal charges are instructional.'
    },
    'TIO2_RUTILE': {
        'label': 'TiO2 crystal fragment (rutile)', 'formula': 'TiO2',
        'prototype': 'rutile, tetragonal P42/mnm', 'formula_units_per_cell': 2,
        'coordination': 'Ti: 6; O: 3',
        'cell': ((1, 0, 0), (0, 1, 0), (0, 0, 0.644)),
        'basis': [
            ('Ti', 4, (0, 0, 0)), ('Ti', 4, (0.5, 0.5, 0.5)),
            ('O', -2, (0.305, 0.305, 0)), ('O', -2, (0.695, 0.695, 0)),
            ('O', -2, (0.805, 0.195, 0.5)), ('O', -2, (0.195, 0.805, 0.5)),
        ],
        'note': 'Idealized rutile topology; Ti(+4)/O(-2) are formal oxidation states, not literal point charges.'
    },
    'CATIO3_PEROVSKITE': {
        'label': 'CaTiO3 crystal fragment (perovskite)', 'formula': 'CaTiO3',
        'prototype': 'cubic perovskite, Pm-3m', 'formula_units_per_cell': 1,
        'coordination': 'Ca: 12; Ti: 6; O: 6 (2 Ti + 4 Ca)',
        'cell': ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        'basis': [
            ('Ca', 2, (0, 0, 0)), ('Ti', 4, (0.5, 0.5, 0.5)),
            ('O', -2, (0.5, 0.5, 0)), ('O', -2, (0.5, 0, 0.5)),
            ('O', -2, (0, 0.5, 0.5)),
        ],
        'note': 'Idealized cubic perovskite; real CaTiO3 distorts at ambient conditions; charges are formal.'
    },
}
IONIC_STRUCTURE_ITEMS = [
    (key, data['label'], f"{data['formula']} — {data['prototype']}; CN {data['coordination']}")
    for key, data in IONIC_STRUCTURE_DATA.items()
]
BOND_ORDER_ITEMS = [
    ('1', "Single", "Single bond"),
    ('2', "Double", "Double bond"),
    ('3', "Triple", "Triple bond"),
]
BUILDER_START_ITEMS = [
    ('CONTINUE', "Continue current model", "Keep the current preset molecule and extend it"),
    ('CLEAR', "Start fresh", "Remove AtomLab-generated model and builder objects first"),
]
BUILDER_MAX_VALENCE = {
    'H': 1, 'He': 0, 'Li': 1, 'Be': 2, 'B': 3, 'C': 4,
    'N': 3, 'O': 2, 'F': 1, 'Ne': 0, 'Na': 1, 'Mg': 2,
    'Al': 3, 'Si': 4, 'P': 5, 'S': 6, 'Cl': 7, 'Ar': 0,
    'K': 1, 'Ca': 2,
}
COVALENT_RADII_A = {
    'H': 0.31, 'Li': 1.28, 'Be': 0.96, 'B': 0.84, 'C': 0.76,
    'N': 0.71, 'O': 0.66, 'F': 0.57, 'Na': 1.66, 'Mg': 1.41,
    'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 1.02,
    'K': 2.03, 'Ca': 1.76,
}


def simple_valence_limit(symbol, formal_charge):
    """Common introductory valence guidelines, not a full bonding theory."""
    charged = {
        'H': {-1: 1, 0: 1, 1: 0},
        'C': {-1: 3, 0: 4, 1: 3},
        'N': {-1: 2, 0: 3, 1: 4},
        'O': {-2: 0, -1: 1, 0: 2, 1: 3},
        'F': {-1: 0, 0: 1, 1: 2},
        'Cl': {-1: 0, 0: 1, 1: 2},
    }
    if symbol in charged and formal_charge in charged[symbol]:
        return charged[symbol][formal_charge]
    return BUILDER_MAX_VALENCE.get(symbol)


# Keep enum strings in a module-level collection. Blender may retain references to
# dynamic EnumProperty item strings; a static list is safer across Blender versions.
ELEMENT_ITEMS = [
    (str(z), f"{symbol} — {name} (Z={z})", f"Neutral {name} atom")
    for z, symbol, name, mass_number in ELEMENTS
]

NUCLEUS_MODE_ITEMS = [
    ('MARKERS', "Nucleon markers", "Show schematic individual proton/neutron markers"),
    ('DENSITY', "Density cloud (approx.)", "Show sampled, illustrative proton/neutron density clouds"),
]


def get_element(z):
    return next(item for item in ELEMENTS if item[0] == z)


def electron_configuration_evidence(z, charge=0):
    """State the evidence tier; ionic configurations use simplified filling rules."""
    if charge:
        return "Ion configuration estimated by outer-shell removal/addition rules."
    if z <= 92:
        return "Reference basis: NIST compilation (H-U)."
    if z <= 103:
        return "NIST ASD/actinide data; some configurations are estimated."
    return "Addon uses Aufbau estimates; heavy-element configurations remain tentative."


def get_subshell_occupancies(z, charge=0):
    """Approximate subshell populations for an atom/ion with charge in e units."""
    if charge > z:
        raise ValueError("Positive charge cannot exceed the atomic number")
    occupancy = {}
    remaining = z
    for n, l, capacity in AUFBAU_ORBITALS:
        count = min(remaining, capacity)
        if count:
            occupancy[(n, l)] = count
            remaining -= count
        if remaining == 0:
            break
    for orbital, change in CONFIGURATION_ADJUSTMENTS.get(z, {}).items():
        occupancy[orbital] = occupancy.get(orbital, 0) + change
    occupancy = {key: count for key, count in occupancy.items() if count > 0}

    if charge > 0:
        electrons_to_remove = charge
        # Ionize the largest principal shell first; this removes ns before (n-1)d,
        # as in the familiar transition-metal cations. Within n, highest l first.
        for orbital in sorted(occupancy, key=lambda item: (item[0], item[1]), reverse=True):
            removed = min(occupancy.get(orbital, 0), electrons_to_remove)
            occupancy[orbital] -= removed
            electrons_to_remove -= removed
            if electrons_to_remove == 0:
                break
        if electrons_to_remove:
            raise ValueError("Charge removes more electrons than the atom contains")
    elif charge < 0:
        electrons_to_add = -charge
        for n, l, capacity in AUFBAU_ORBITALS:
            current = occupancy.get((n, l), 0)
            added = min(capacity - current, electrons_to_add)
            if added > 0:
                occupancy[(n, l)] = current + added
                electrons_to_add -= added
            if electrons_to_add == 0:
                break
        if electrons_to_add:
            raise ValueError("Could not place all additional electrons in configured subshells")

    return {key: count for key, count in occupancy.items() if count > 0}


def format_electron_configuration(z, charge=0):
    """Return a compact approximate Aufbau-order configuration for an atom or ion."""
    occupancy = get_subshell_occupancies(z, charge)
    labels = 'spdf'
    return " ".join(
        f"{n}{labels[l]}{occupancy[(n, l)]}"
        for n, l, _capacity in AUFBAU_ORBITALS
        if occupancy.get((n, l), 0) > 0
    )


def get_shell_occupancies(z, charge=0):
    """Return shell populations summed from the approximate configuration."""
    occupancy = {}
    for (n, _l), count in get_subshell_occupancies(z, charge).items():
        occupancy[n] = occupancy.get(n, 0) + count
    return sorted(occupancy.items())


def slater_effective_charge(z, subshells, n, l, target_electron_included=True):
    """Simple Slater-style Zeff estimate; suitable only as a teaching approximation."""
    # For a one-electron atom, hydrogen-like orbitals feel the full nuclear charge.
    if z == 1:
        return 1.0
    # In a simple He 1s -> n=2 excited-orbital picture, one 1s electron remains.
    if z == 2 and n == 2 and not target_electron_included:
        return 1.15
    if n == 1 and l == 0:
        shielding = 0.30 * max(0, subshells.get((1, 0), 0) - int(target_electron_included))
    elif l >= 2:
        same_group = subshells.get((n, l), 0) - int(target_electron_included)
        shielding = 0.35 * max(0, same_group)
        # For d/f orbitals, groups to their left (lower n or lower l at the same n)
        # shield approximately one-for-one in the Slater rules.
        for (other_n, other_l), count in subshells.items():
            if (other_n, other_l) == (n, l):
                continue
            if other_n < n or (other_n == n and other_l < l):
                shielding += count
    else:
        same_group = subshells.get((n, 0), 0) + subshells.get((n, 1), 0)
        shielding = 0.35 * max(0, same_group - int(target_electron_included))
        for (other_n, other_l), count in subshells.items():
            if other_n == n:
                continue
            if other_n == n - 1:
                shielding += 0.85 * count
            elif other_n < n - 1:
                shielding += count
    return max(0.5, z - shielding)


def associated_laguerre(order, alpha, x):
    if order == 0:
        return 1.0
    previous_previous = 1.0
    previous = alpha + 1.0 - x
    if order == 1:
        return previous
    for k in range(2, order + 1):
        current = ((2 * k - 1 + alpha - x) * previous
                   - (k - 1 + alpha) * previous_previous) / k
        previous_previous, previous = previous, current
    return previous


def orbital_angular_weight(l, orientation, x, y, z):
    """Peak-scaled real spherical-harmonic angular probability factors for s/p/d/f."""
    if l == 0:
        return 1.0
    if l == 1:
        component = {'X': x, 'Y': y, 'Z': z}[orientation]
        return component * component
    if l == 2:
        weights = {
            'DXY': 4.0 * x*x * y*y,
            'DXZ': 4.0 * x*x * z*z,
            'DYZ': 4.0 * y*y * z*z,
            'DX2Y2': (x*x - y*y) ** 2,
            'DZ2': ((3.0 * z*z - 1.0) / 2.0) ** 2,
        }
        return min(1.0, weights[orientation])
    if l == 3:
        weights = {
            'FZ3': (z * (5.0 * z*z - 3.0) / 2.0) ** 2,
            'FXZ2': (x * (5.0 * z*z - 1.0) / 4.0) ** 2,
            'FYZ2': (y * (5.0 * z*z - 1.0) / 4.0) ** 2,
            'FZXY': (z * (x*x - y*y)) ** 2,
            'FXYZ': (3.0 * math.sqrt(3.0) * x * y * z) ** 2,
            'FXX23Y2': (x * (x*x - 3.0 * y*y) / 3.0) ** 2,
            'FY3X2Y2': (y * (3.0 * x*x - y*y) / 3.0) ** 2,
        }
        return min(1.0, weights[orientation])
    raise ValueError(f"Unsupported angular momentum l={l}")


def sample_hydrogenic_orbital(n, l, orientation, zeff, count, seed):
    """Sample a hydrogen-like probability density; multi-electron Zeff is approximate."""
    rng = random.Random(seed)
    radial_order = n - l - 1
    alpha = 2 * l + 1
    r_max = max(12.0, 14.0 * n * n / zeff)
    steps = 1800
    dr = r_max / steps
    radii = [i * dr for i in range(steps + 1)]
    weights = []
    for r in radii:
        rho = 2.0 * zeff * r / n
        laguerre = associated_laguerre(radial_order, alpha, rho)
        weight = (r ** (2 * l + 2)) * math.exp(-2.0 * zeff * r / n) * laguerre ** 2
        weights.append(weight)

    cdf = [0.0]
    for i in range(1, len(weights)):
        cdf.append(cdf[-1] + 0.5 * (weights[i - 1] + weights[i]) * dr)
    total = cdf[-1]
    if total <= 0.0:
        raise ValueError("Could not calculate orbital radial distribution")

    points = []
    while len(points) < count:
        target = rng.random() * total
        low, high = 0, len(cdf) - 1
        while low < high:
            middle = (low + high) // 2
            if cdf[middle] < target:
                low = middle + 1
            else:
                high = middle
        index = max(1, low)
        before, after = cdf[index - 1], cdf[index]
        fraction = 0.0 if after == before else (target - before) / (after - before)
        r = (radii[index - 1] + fraction * dr) * 0.48  # visual Bohr-radius scale

        while True:
            z_dir = rng.uniform(-1.0, 1.0)
            phi = rng.uniform(0.0, 2.0 * math.pi)
            planar = math.sqrt(max(0.0, 1.0 - z_dir * z_dir))
            x_dir = planar * math.cos(phi)
            y_dir = planar * math.sin(phi)
            if l == 0 or rng.random() <= orbital_angular_weight(
                l, orientation, x_dir, y_dir, z_dir
            ):
                break
        points.append((r * x_dir, r * y_dir, r * z_dir))
    return points


def get_or_create_emission_material(name, color, strength):
    mat = bpy.data.materials.get(name)
    if mat is None:
        mat = bpy.data.materials.new(name=name)
    # Blender 5+ creates node-enabled materials by default; setting use_nodes
    # is deprecated there (and scheduled for removal in Blender 6.0).
    if bpy.app.version < (5, 0, 0):
        mat.use_nodes = True
    mat.diffuse_color = (*color, 1.0)
    nodes = mat.node_tree.nodes
    links = mat.node_tree.links
    nodes.clear()
    emission = nodes.new("ShaderNodeEmission")
    emission.inputs["Color"].default_value = (*color, 1.0)
    emission.inputs["Strength"].default_value = strength
    output = nodes.new("ShaderNodeOutputMaterial")
    links.new(emission.outputs["Emission"], output.inputs["Surface"])
    return mat


def put_in_collection(obj, collection):
    for old_collection in tuple(obj.users_collection):
        old_collection.objects.unlink(obj)
    collection.objects.link(obj)


def remove_previous_atomlab_model():
    # Also remove older suffixed copies (e.g. AtomLab_Generated.001) left by
    # previous prototype runs, while preserving the separate Builder collection.
    generated = [
        collection for collection in bpy.data.collections
        if collection.name == "AtomLab_Generated"
        or collection.name.startswith("AtomLab_Generated.")
        or collection.get("atomlab_managed_model", False)
    ]
    for collection in generated:
        for obj in tuple(collection.objects):
            bpy.data.objects.remove(obj, do_unlink=True)
        bpy.data.collections.remove(collection, do_unlink=True)


def create_model_collection(scene):
    remove_previous_atomlab_model()
    collection = bpy.data.collections.new("AtomLab_Generated")
    collection["atomlab_managed_model"] = True
    scene.collection.children.link(collection)
    return collection


def add_sphere(name, radius, location, collection, material):
    bpy.ops.mesh.primitive_uv_sphere_add(
        segments=32, ring_count=16, radius=radius, location=location
    )
    obj = bpy.context.object
    obj.name = name
    for polygon in obj.data.polygons:
        polygon.use_smooth = True
    if material:
        obj.data.materials.append(material)
    put_in_collection(obj, collection)
    return obj


def nucleon_positions(total, extent):
    """Deterministic positions inside a schematic nucleus."""
    if total == 1:
        return [(0.0, 0.0, 0.0)]
    golden_angle = math.pi * (3.0 - math.sqrt(5.0))
    positions = []
    for i in range(total):
        z_unit = 1.0 - 2.0 * (i + 0.5) / total
        ring = math.sqrt(max(0.0, 1.0 - z_unit * z_unit))
        angle = i * golden_angle
        radial = extent * ((i + 0.5) / total) ** (1.0 / 3.0)
        positions.append((
            radial * ring * math.cos(angle),
            radial * ring * math.sin(angle),
            radial * z_unit,
        ))
    return positions


def add_density_cloud(name, nucleon_count, extent, collection, material, seed, isotope):
    """Create many tiny samples of an illustrative, isotropic nuclear density."""
    if nucleon_count <= 0:
        return None

    rng = random.Random(seed)
    sample_count = max(24, nucleon_count * 18)
    sigma = extent * 0.44
    dot_radius = 0.010
    vertices = []
    faces = []
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    ico = [
        (-1, phi, 0), (1, phi, 0), (-1, -phi, 0), (1, -phi, 0),
        (0, -1, phi), (0, 1, phi), (0, -1, -phi), (0, 1, -phi),
        (phi, 0, -1), (phi, 0, 1), (-phi, 0, -1), (-phi, 0, 1),
    ]
    ico = [tuple(c / math.sqrt(x*x + y*y + z*z) for c in (x, y, z))
           for x, y, z in ico]
    ico_faces = (
        (0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
        (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
        (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
        (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1),
    )

    for _ in range(sample_count):
        while True:
            x = rng.gauss(0.0, sigma)
            y = rng.gauss(0.0, sigma)
            z = rng.gauss(0.0, sigma)
            if x*x + y*y + z*z <= extent*extent:
                break
        base = len(vertices)
        vertices.extend((x + dot_radius * vx, y + dot_radius * vy, z + dot_radius * vz)
                        for vx, vy, vz in ico)
        faces.extend(tuple(base + v for v in face) for face in ico_faces)

    mesh = bpy.data.meshes.new(name + "_Mesh")
    mesh.from_pydata(vertices, [], faces)
    mesh.update()
    for polygon in mesh.polygons:
        polygon.use_smooth = True
    obj = bpy.data.objects.new(name, mesh)
    collection.objects.link(obj)
    mesh.materials.append(material)
    obj["sample_count"] = sample_count
    obj["isotope"] = isotope
    obj["educational_note"] = (
        "Illustrative isotropic density samples; not measured nucleon positions "
        "or an isotope-specific nuclear wavefunction."
    )
    return obj


def add_orbital_point_mesh(name, points, collection, material, orbital, occupancy, zeff):
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    base_vertices = [
        (-1, phi, 0), (1, phi, 0), (-1, -phi, 0), (1, -phi, 0),
        (0, -1, phi), (0, 1, phi), (0, -1, -phi), (0, 1, -phi),
        (phi, 0, -1), (phi, 0, 1), (-phi, 0, -1), (-phi, 0, 1),
    ]
    base_vertices = [
        tuple(c / math.sqrt(x*x + y*y + z*z) for c in (x, y, z))
        for x, y, z in base_vertices
    ]
    base_faces = (
        (0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
        (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
        (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
        (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1),
    )
    dot_radius = 0.011
    vertices, faces = [], []
    for x, y, z in points:
        offset = len(vertices)
        vertices.extend((x + dot_radius * vx, y + dot_radius * vy, z + dot_radius * vz)
                        for vx, vy, vz in base_vertices)
        faces.extend(tuple(offset + i for i in face) for face in base_faces)

    mesh = bpy.data.meshes.new(name + "_Mesh")
    mesh.from_pydata(vertices, [], faces)
    mesh.update()
    for polygon in mesh.polygons:
        polygon.use_smooth = True
    obj = bpy.data.objects.new(name, mesh)
    collection.objects.link(obj)
    mesh.materials.append(material)
    obj["orbital"] = orbital
    obj["subshell_electron_count"] = occupancy
    obj["effective_nuclear_charge_approx"] = zeff
    obj["educational_note"] = (
        "Single-orbital hydrogen-like probability samples with approximate Zeff. "
        "The recorded electron count belongs to the subshell, not this one real "
        "orientation. d/f orientations are representative real-harmonic basis shapes."
    )
    return obj


def animate_electron_schematic(obj, radius, initial_angle, period_frames):
    """Animate an electron directly, avoiding visible parent/child guide lines."""
    for axis, function in ((0, "cos"), (1, "sin")):
        fcurve = obj.driver_add("location", axis)
        driver = fcurve.driver
        driver.type = 'SCRIPTED'
        driver.expression = (
            f"{radius:.9f} * {function}({initial_angle:.9f} + "
            f"6.283185307179586 * (frame - 1) / {period_frames:.9f})"
        )


def add_orbit(name, radius, collection, material):
    bpy.ops.mesh.primitive_torus_add(
        major_radius=radius,
        minor_radius=0.012,
        major_segments=96,
        minor_segments=8,
        location=(0.0, 0.0, 0.0),
    )
    obj = bpy.context.object
    obj.name = name
    for polygon in obj.data.polygons:
        polygon.use_smooth = True
    obj.data.materials.append(material)
    put_in_collection(obj, collection)
    return obj


def ionic_structure_sites(structure_id, repeats=(2, 2, 2)):
    """Expand fractional basis sites into a centered finite crystal fragment."""
    if structure_id not in IONIC_STRUCTURE_DATA:
        raise ValueError(f"Unsupported ionic structure: {structure_id}")
    data = IONIC_STRUCTURE_DATA[structure_id]
    if len(repeats) != 3 or any(int(count) < 1 for count in repeats):
        raise ValueError("Crystal repeats must be three positive integers")
    repeats = tuple(int(count) for count in repeats)
    cell = data['cell']
    raw_sites = []
    for i in range(repeats[0]):
        for j in range(repeats[1]):
            for k in range(repeats[2]):
                cell_offset = (i, j, k)
                for element, charge, fractional in data['basis']:
                    coordinates = tuple(
                        sum((cell_offset[axis] + fractional[axis]) * cell[axis][component]
                            for axis in range(3))
                        for component in range(3)
                    )
                    raw_sites.append((element, charge, coordinates))
    center = tuple(
        (min(site[2][axis] for site in raw_sites)
         + max(site[2][axis] for site in raw_sites)) / 2.0
        for axis in range(3)
    )
    return [
        (element, charge, tuple(point[axis] - center[axis] for axis in range(3)))
        for element, charge, point in raw_sites
    ]


def create_ionic_structure(context, structure_id):
    if structure_id not in IONIC_STRUCTURE_DATA:
        raise ValueError(f"Unsupported ionic structure: {structure_id}")
    data = IONIC_STRUCTURE_DATA[structure_id]
    repeats = (2, 2, 2)
    sites = ionic_structure_sites(structure_id, repeats)
    total_charge = sum(charge for _element, charge, _position in sites)
    if total_charge != 0:
        raise ValueError(f"Ionic fragment is not charge-neutral (net charge {total_charge})")

    scene = context.scene
    collection = create_model_collection(scene)
    formula_units = data['formula_units_per_cell'] * math.prod(repeats)
    collection["structure_type"] = "Finite ionic crystal fragment"
    collection["formula_units"] = formula_units
    collection["crystal_structure"] = data['prototype']

    scale = 2.7  # schematic Blender units per normalized cell edge
    colors = {
        "Na": (1.0, 0.58, 0.08), "Cl": (0.10, 0.78, 0.36),
        "Cs": (0.66, 0.36, 0.95), "Zn": (0.12, 0.48, 1.0),
        "S": (1.0, 0.82, 0.12), "Ca": (1.0, 0.40, 0.08),
        "F": (0.05, 0.88, 0.82), "Li": (0.95, 0.35, 0.62),
        "O": (1.0, 0.16, 0.12), "Ti": (0.38, 0.56, 0.82),
        "Ni": (0.62, 0.68, 0.76), "As": (0.72, 0.34, 0.86),
        "Cd": (0.35, 0.76, 0.68), "I": (0.56, 0.29, 0.76),
        "Mg": (0.72, 0.74, 0.79), "Al": (0.92, 0.65, 0.26),
    }
    radii = {
        "Na": 0.46, "Cl": 0.66, "Cs": 0.62, "Zn": 0.43, "S": 0.58,
        "Ca": 0.46, "F": 0.42, "Li": 0.34, "O": 0.45, "Ti": 0.48,
        "Ni": 0.47, "As": 0.55, "Cd": 0.52, "I": 0.60,
        "Mg": 0.25, "Al": 0.19,
    }
    materials = {
        element: get_or_create_emission_material(
            f"AtomLab_Ionic_{element}_Mat", color, 1.1
        )
        for element, color in colors.items()
        if any(site[0] == element for site in sites)
    }
    counts = collections.Counter()
    for element, charge, position in sites:
        counts[element] += 1
        symbol_charge = f"{charge:+d}"
        location = tuple(coord * scale for coord in position)
        ion = add_sphere(
            f"{structure_id}_{element}_{counts[element]:03d}",
            radii[element],
            location,
            collection,
            materials[element],
        )
        element_row = next(item for item in ELEMENTS if item[1] == element)
        ion["element_symbol"] = element
        ion["atomic_number"] = element_row[0]
        ion["element_name"] = element_row[2]
        ion["ionic_species"] = f"{element}{symbol_charge}"
        ion["formal_charge"] = charge
        ion["crystal_structure"] = data['prototype']
        ion["bulk_coordination_summary"] = data['coordination']
        ion["finite_fragment_surface"] = True
        ion["is_discrete_molecule"] = False
        ion["educational_note"] = (
            f"{data['formula']} is a finite {data['prototype']} crystal fragment, "
            "not separate molecular units. Formal integer charges are teaching labels; "
            "surface coordination is reduced, and no covalent rods are drawn. "
            + data['note']
        )

    scene["atomlab_model"] = "Ionic crystal fragment"
    scene["atomlab_ionic_structure_id"] = structure_id
    scene["atomlab_formula"] = data['formula']
    scene["atomlab_formula_unit_count"] = formula_units
    scene["atomlab_ionic_prototype"] = data['prototype']
    scene["atomlab_ionic_note"] = (
        f"Finite {data['prototype']} fragment, not discrete {data['formula']} molecules. "
        f"Bulk coordination: {data['coordination']}. {data['note']}"
    )
    scene["atomlab_builder_net_charge"] = builder_net_charge(scene)
    return {
        "formula": data['formula'],
        "formula_units": formula_units,
        "ion_counts": dict(counts),
        "structure": data['prototype'],
        "coordination": data['coordination'],
    }


def create_molecule_model(context, molecule_id):
    data = MOLECULE_DATA[molecule_id]
    scene = context.scene
    collection = create_model_collection(scene)
    scale = 1.4  # Blender units per Angstrom, chosen for comfortable viewport framing.

    colors = {
        'H': (0.82, 0.90, 1.00),
        'C': (0.13, 0.16, 0.21),
        'N': (0.10, 0.34, 1.00),
        'O': (1.00, 0.12, 0.10),
    }
    atom_materials = {
        element: get_or_create_emission_material(
            f"AtomLab_Molecule_{element}_Mat", color, 0.8
        )
        for element, color in colors.items()
    }
    bond_material = get_or_create_emission_material(
        "AtomLab_Molecule_Bond_Mat", (0.58, 0.68, 0.82), 0.55
    )
    atom_radii = {'H': 0.18, 'C': 0.29, 'N': 0.28, 'O': 0.28}
    atom_objects = []

    for index, (element, coordinates) in enumerate(data['atoms']):
        location = Vector(coordinates) * scale
        atom = add_sphere(
            f"{molecule_id}_Atom_{index + 1}_{element}",
            atom_radii[element],
            tuple(location),
            collection,
            atom_materials[element],
        )
        atom["element"] = element
        atom["element_symbol"] = element
        atom["atomic_number"] = next(item[0] for item in ELEMENTS if item[1] == element)
        atom["element_name"] = next(item[2] for item in ELEMENTS if item[1] == element)
        atom["molecule_formula"] = data['formula']
        atom["atomlab_builder_atom"] = True
        atom["formal_charge"] = 0
        atom["atomlab_bond_order_sum"] = sum(
            order for a, b, order in data['bonds'] if index in (a, b)
        )
        atom_objects.append(atom)

    for bond_index, (first, second, order) in enumerate(data['bonds'], start=1):
        start = Vector(data['atoms'][first][1]) * scale
        end = Vector(data['atoms'][second][1]) * scale
        direction = (end - start).normalized()
        length = (end - start).length
        reference = Vector((0.0, 0.0, 1.0))
        if abs(direction.dot(reference)) > 0.90:
            reference = Vector((0.0, 1.0, 0.0))
        perpendicular = direction.cross(reference).normalized()
        if order == 1:
            offsets = (0.0,)
        elif order == 2:
            offsets = (-0.045, 0.045)
        else:
            offsets = (-0.07, 0.0, 0.07)

        for strand, offset in enumerate(offsets, start=1):
            bpy.ops.mesh.primitive_cylinder_add(
                vertices=16,
                radius=0.032,
                depth=length,
                location=tuple((start + end) * 0.5 + perpendicular * offset),
            )
            rod = bpy.context.object
            rod.name = f"{molecule_id}_Bond_{bond_index}_{strand}"
            rod.rotation_euler = direction.to_track_quat('Z', 'Y').to_euler()
            rod.data.materials.append(bond_material)
            for polygon in rod.data.polygons:
                polygon.use_smooth = True
            put_in_collection(rod, collection)
            rod["bond_order"] = order
            rod["bond_length_angstrom"] = data['bond_length_A']
            rod["bond_base_length"] = length
            rod["bond_offset"] = offset
            rod["atomlab_builder_bond"] = True
            rod["atom_a_name"] = atom_objects[first].name
            rod["atom_b_name"] = atom_objects[second].name

    scene["atomlab_model"] = "Molecule library (ball-and-stick)"
    scene["atomlab_molecule_name"] = data['label']
    scene["atomlab_molecule_id"] = molecule_id
    scene["atomlab_formula"] = data['formula']
    scene["atomlab_geometry"] = data['geometry']
    scene["atomlab_bond_lengths_angstrom"] = data['bond_length_A']
    scene["atomlab_builder_net_charge"] = builder_net_charge(scene)
    return data


def get_or_create_builder_collection(scene):
    collection = bpy.data.collections.get("AtomLab_Builder")
    if collection is None:
        collection = bpy.data.collections.new("AtomLab_Builder")
    if collection.name not in {child.name for child in scene.collection.children}:
        scene.collection.children.link(collection)
    return collection


def tag_current_preset_for_builder(scene):
    """Tag a visible preset molecule from an earlier addon run as builder-compatible."""
    generated = bpy.data.collections.get("AtomLab_Generated")
    if generated is None or scene.get("atomlab_model") != "Molecule library (ball-and-stick)":
        return
    molecule_id = scene.get("atomlab_molecule_id")
    if molecule_id not in MOLECULE_DATA:
        formula = scene.get("atomlab_formula")
        molecule_id = next(
            (key for key, data in MOLECULE_DATA.items() if data['formula'] == formula),
            None,
        )
    if molecule_id not in MOLECULE_DATA:
        return

    data = MOLECULE_DATA[molecule_id]
    atom_objects = {}
    for index, (element, _coordinates) in enumerate(data['atoms']):
        expected_name = f"{molecule_id}_Atom_{index + 1}_{element}"
        atom = next(
            (obj for obj in generated.objects
             if obj.name == expected_name or obj.name.startswith(expected_name + ".")),
            None,
        )
        if atom is None:
            continue
        atomic_number, _symbol, name, _mass = next(item for item in ELEMENTS if item[1] == element)
        atom["atomlab_builder_atom"] = True
        atom["atomic_number"] = atomic_number
        atom["element_symbol"] = element
        atom["element_name"] = name
        atom["formal_charge"] = 0
        atom["atomlab_bond_order_sum"] = sum(
            order for a, b, order in data['bonds'] if index in (a, b)
        )
        atom_objects[index] = atom

    for bond_index, (first, second, order) in enumerate(data['bonds'], start=1):
        if first not in atom_objects or second not in atom_objects:
            continue
        prefix = f"{molecule_id}_Bond_{bond_index}_"
        for rod in generated.objects:
            if rod.name.startswith(prefix):
                rod["atomlab_builder_bond"] = True
                rod["atom_a_name"] = atom_objects[first].name
                rod["atom_b_name"] = atom_objects[second].name
                rod["bond_order"] = order


def prepare_builder(context, start_mode):
    if start_mode == 'CLEAR':
        remove_previous_atomlab_model()
        clear_builder_collection(context.scene)
    else:
        tag_current_preset_for_builder(context.scene)
    collection = get_or_create_builder_collection(context.scene)
    collection["atomlab_builder_session_started"] = True
    context.scene["atomlab_builder_net_charge"] = builder_net_charge(context.scene)
    return collection


def add_builder_atom(context, z, formal_charge=0):
    collection = get_or_create_builder_collection(context.scene)
    atomic_number, symbol, name, mass_number = get_element(z)
    existing = [obj for obj in collection.objects if obj.get("atomlab_builder_atom", False)]
    index = len(existing) + 1
    # Put new atoms to the right of the current fragment, then let the user refine with G.
    all_atoms = [obj for obj in context.scene.objects if obj.get("atomlab_builder_atom", False)]
    cursor = context.scene.cursor.location.copy()
    if all_atoms:
        rightmost = max(obj.matrix_world.translation.x for obj in all_atoms)
        location = Vector((rightmost + 1.25, cursor.y, cursor.z))
    else:
        location = cursor
    colors = {
        'H': (0.82, 0.90, 1.00), 'C': (0.13, 0.16, 0.21),
        'N': (0.10, 0.34, 1.00), 'O': (1.00, 0.12, 0.10),
        'S': (1.00, 0.72, 0.10), 'P': (1.00, 0.45, 0.10),
        'F': (0.35, 0.90, 0.30), 'Cl': (0.25, 0.85, 0.30),
    }
    color = colors.get(symbol, (0.95, 0.55, 0.12))
    material = get_or_create_emission_material(
        f"AtomLab_Builder_{symbol}_Mat", color, 0.85
    )
    radius = {'H': 0.18, 'C': 0.29, 'N': 0.28, 'O': 0.28}.get(symbol, 0.25)
    atom = add_sphere(
        f"Builder_{symbol}_{index:02d}", radius, tuple(location), collection, material
    )
    atom["atomlab_builder_atom"] = True
    atom["atomic_number"] = atomic_number
    atom["element_symbol"] = symbol
    atom["element_name"] = name
    atom["formal_charge"] = int(formal_charge)
    atom["atomlab_bond_order_sum"] = 0
    atom["educational_note"] = "Builder atom marker; move before adding bonds."
    context.scene["atomlab_builder_net_charge"] = builder_net_charge(context.scene)
    return atom


def builder_net_charge(scene):
    return sum(
        int(obj.get("formal_charge", 0))
        for obj in scene.objects
        if obj.get("atomlab_builder_atom", False)
    )


def set_selected_atom_charge(context, formal_charge):
    selected = [obj for obj in context.selected_objects
                if obj.get("atomlab_builder_atom", False)]
    if len(selected) != 1:
        raise ValueError("Select exactly one AtomLab atom before setting its charge.")
    atom = selected[0]
    charge = int(formal_charge)
    atom["formal_charge"] = charge
    valence = int(atom.get("atomlab_bond_order_sum", 0))
    limit = simple_valence_limit(atom["element_symbol"], charge)
    warnings = []
    if limit is not None and valence > limit:
        warnings.append(
            f"current bond-order sum {valence} exceeds the simple guideline {limit} "
            f"for charge {charge:+d}"
        )
    total_charge = builder_net_charge(context.scene)
    context.scene["atomlab_builder_net_charge"] = total_charge
    return atom, total_charge, warnings


def empty_custom_molecule_library():
    return {'schema_version': CUSTOM_MOLECULE_SCHEMA_VERSION, 'molecules': []}


def custom_library_user_path():
    directory = bpy.utils.user_resource('CONFIG', path='atomlab', create=True)
    if not directory:
        directory = os.path.join(bpy.utils.resource_path('USER'), 'config', 'atomlab')
        os.makedirs(directory, exist_ok=True)
    return os.path.join(directory, 'custom_molecules.json')


def read_project_custom_library():
    text = bpy.data.texts.get(CUSTOM_MOLECULE_TEXT_NAME)
    if text is None or not text.as_string().strip():
        return empty_custom_molecule_library()
    return validate_custom_library_payload(json.loads(text.as_string()))


def write_project_custom_library(payload):
    payload = validate_custom_library_payload(payload)
    text = bpy.data.texts.get(CUSTOM_MOLECULE_TEXT_NAME)
    if text is None:
        text = bpy.data.texts.new(CUSTOM_MOLECULE_TEXT_NAME)
    text.clear()
    text.write(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True))
    text.use_fake_user = True


def read_global_custom_library():
    path = custom_library_user_path()
    if not os.path.exists(path):
        return empty_custom_molecule_library()
    with open(path, 'r', encoding='utf-8') as handle:
        return validate_custom_library_payload(json.load(handle))


def write_global_custom_library(payload):
    payload = validate_custom_library_payload(payload)
    path = custom_library_user_path()
    directory = os.path.dirname(path)
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode='w', encoding='utf-8', dir=directory, prefix='.atomlab-',
            suffix='.json.tmp', delete=False,
        ) as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.write('\n')
            temporary_path = handle.name
        os.replace(temporary_path, path)
    finally:
        if temporary_path and os.path.exists(temporary_path):
            os.remove(temporary_path)
    return path


def refresh_custom_molecule_catalog():
    """Merge project-embedded and machine-wide JSON records for the UI."""
    sources = {}
    records = {}
    project = read_project_custom_library()
    shared = read_global_custom_library()
    for label, payload in (("Project", project), ("Shared", shared)):
        for record in payload['molecules']:
            key = 'CUSTOM_' + uuid.UUID(record['id']).hex.upper()
            if key in records:
                if records[key] != record and label == "Project":
                    records[key] = record
                sources[key].add(label)
            else:
                records[key] = record
                sources[key] = {label}
    CUSTOM_MOLECULES_BY_ID.clear()
    CUSTOM_MOLECULES_BY_ID.update(records)
    items = [(CUSTOM_MOLECULE_NONE, "No saved molecules", "Save a molecule from the Builder first")]
    for key, record in sorted(records.items(), key=lambda pair: pair[1]['name'].casefold()):
        source_label = " + ".join(sorted(sources[key]))
        description = f"{format_library_display_text(record['formula'])} | charge {record['net_charge']:+d} | {source_label}"
        items.append((key, format_library_display_text(record['name']), description[:240]))
    CUSTOM_MOLECULE_ITEMS[:] = items
    return len(records)


def custom_molecule_enum_items(_scene, _context):
    return CUSTOM_MOLECULE_ITEMS


def custom_library_payload_for_export():
    project = read_project_custom_library()
    shared = read_global_custom_library()
    combined = {}
    for payload in (shared, project):
        for record in payload['molecules']:
            combined[record['id']] = record
    return {
        'schema_version': CUSTOM_MOLECULE_SCHEMA_VERSION,
        'molecules': sorted(combined.values(), key=lambda record: record['name'].casefold()),
    }


def _write_custom_libraries_transactionally(project_payload, shared_payload):
    text = bpy.data.texts.get(CUSTOM_MOLECULE_TEXT_NAME)
    old_project_text = text.as_string() if text is not None else None
    global_path = custom_library_user_path()
    old_shared_bytes = None
    if os.path.isfile(global_path):
        with open(global_path, 'rb') as handle:
            old_shared_bytes = handle.read()
    try:
        write_project_custom_library(project_payload)
        write_global_custom_library(shared_payload)
    except Exception:
        # Best-effort compensating rollback so a failed second write does not leave
        # the project and shared catalogs silently out of sync.
        if old_project_text is None:
            stale = bpy.data.texts.get(CUSTOM_MOLECULE_TEXT_NAME)
            if stale is not None:
                bpy.data.texts.remove(stale)
        else:
            restored_text = bpy.data.texts.get(CUSTOM_MOLECULE_TEXT_NAME)
            if restored_text is None:
                restored_text = bpy.data.texts.new(CUSTOM_MOLECULE_TEXT_NAME)
            restored_text.clear()
            restored_text.write(old_project_text)
            restored_text.use_fake_user = True
        if old_shared_bytes is None:
            try:
                os.remove(global_path)
            except FileNotFoundError:
                pass
        else:
            directory = os.path.dirname(global_path)
            temporary_path = None
            try:
                with tempfile.NamedTemporaryFile(
                    mode='wb', dir=directory, prefix='.atomlab-rollback-',
                    suffix='.json.tmp', delete=False,
                ) as handle:
                    handle.write(old_shared_bytes)
                    handle.flush()
                    os.fsync(handle.fileno())
                    temporary_path = handle.name
                os.replace(temporary_path, global_path)
            finally:
                if temporary_path and os.path.exists(temporary_path):
                    os.remove(temporary_path)
        raise
    return global_path


def save_custom_library_record_to_both(context, record, replace_conflicts=False):
    validate_custom_molecule_record(record)
    project = read_project_custom_library()
    shared = read_global_custom_library()
    project_records = merge_custom_library_records(
        project['molecules'], [record], replace_conflicts
    )
    shared_records = merge_custom_library_records(
        shared['molecules'], [record], replace_conflicts
    )
    new_project = {'schema_version': CUSTOM_MOLECULE_SCHEMA_VERSION, 'molecules': project_records}
    new_shared = {'schema_version': CUSTOM_MOLECULE_SCHEMA_VERSION, 'molecules': shared_records}
    global_path = _write_custom_libraries_transactionally(new_project, new_shared)
    refresh_custom_molecule_catalog()
    return global_path


def import_custom_library_to_both(payload, replace_conflicts=False):
    payload = validate_custom_library_payload(payload)
    project = read_project_custom_library()
    shared = read_global_custom_library()
    project_records = merge_custom_library_records(
        project['molecules'], payload['molecules'], replace_conflicts
    )
    shared_records = merge_custom_library_records(
        shared['molecules'], payload['molecules'], replace_conflicts
    )
    new_project = {'schema_version': CUSTOM_MOLECULE_SCHEMA_VERSION, 'molecules': project_records}
    new_shared = {'schema_version': CUSTOM_MOLECULE_SCHEMA_VERSION, 'molecules': shared_records}
    global_path = _write_custom_libraries_transactionally(new_project, new_shared)
    refresh_custom_molecule_catalog()
    return global_path, len(payload['molecules'])


def custom_record_from_builder(context, name, note=''):
    atoms = sorted(
        (obj for obj in context.scene.objects if obj.get('atomlab_builder_atom', False)),
        key=lambda obj: (obj.name.casefold(), obj.name),
    )
    if not atoms:
        raise ValueError("No AtomLab Builder atoms found; prepare or build a molecule first.")
    atom_index = {atom.name: index for index, atom in enumerate(atoms)}
    world_positions = [atom.matrix_world.translation.copy() for atom in atoms]
    center = sum(world_positions, Vector((0.0, 0.0, 0.0))) / len(world_positions)
    atom_records = []
    for atom, world_position in zip(atoms, world_positions):
        symbol = atom.get('element_symbol')
        if not symbol:
            raise ValueError(f"Builder atom {atom.name} has no element symbol.")
        position = world_position - center
        atom_records.append({
            'element': str(symbol),
            'charge': int(atom.get('formal_charge', 0)),
            'position': [round(float(value), 7) for value in position],
        })
    bond_map = {}
    for obj in context.scene.objects:
        if not obj.get('atomlab_builder_bond', False):
            continue
        first_name, second_name = obj.get('atom_a_name'), obj.get('atom_b_name')
        if first_name not in atom_index or second_name not in atom_index:
            continue
        first, second = atom_index[first_name], atom_index[second_name]
        pair = tuple(sorted((first, second)))
        order = int(obj.get('bond_order', 1))
        previous = bond_map.get(pair)
        if previous is not None and previous != order:
            raise ValueError("Bond rods disagree about bond order; repair the Builder bonds before saving.")
        bond_map[pair] = order
    bond_records = [
        {'a': first, 'b': second, 'order': order}
        for (first, second), order in sorted(bond_map.items())
    ]
    symbols = [atom['element'] for atom in atom_records]
    record = {
        'schema_version': CUSTOM_MOLECULE_SCHEMA_VERSION,
        'id': uuid.uuid4().hex,
        'name': name.strip(),
        'formula': custom_formula_from_symbols(symbols),
        'coordinate_unit': 'blender_unit',
        'atoms': atom_records,
        'bonds': bond_records,
        'net_charge': sum(atom['charge'] for atom in atom_records),
        'note': str(note).strip(),
    }
    return validate_custom_molecule_record(record)


def create_saved_custom_molecule(context, record):
    validate_custom_molecule_record(record)
    scene = context.scene
    collection = create_model_collection(scene)
    collection.name = "AtomLab_Custom_" + uuid.UUID(record['id']).hex[:8]
    collection['atomlab_custom_molecule_id'] = record['id']
    collection['molecule_name'] = record['name']
    collection['molecule_formula'] = record['formula']
    collection['net_charge'] = record['net_charge']

    colors = {
        'H': (0.82, 0.90, 1.00), 'C': (0.13, 0.16, 0.21),
        'N': (0.10, 0.34, 1.00), 'O': (1.00, 0.12, 0.10),
        'S': (1.00, 0.72, 0.10), 'P': (1.00, 0.45, 0.10),
        'F': (0.35, 0.90, 0.30), 'Cl': (0.25, 0.85, 0.30),
    }
    radii = {'H': 0.18, 'C': 0.29, 'N': 0.28, 'O': 0.28}
    elements = {item[1]: item for item in ELEMENTS}
    materials = {}
    atom_objects = []
    for symbol in {atom['element'] for atom in record['atoms']}:
        materials[symbol] = get_or_create_emission_material(
            f"AtomLab_Custom_{symbol}_Mat", colors.get(symbol, (0.95, 0.55, 0.12)), 0.8
        )
    bond_material = get_or_create_emission_material(
        "AtomLab_Custom_Bond_Mat", (0.58, 0.68, 0.82), 0.55
    )
    for index, atom_record in enumerate(record['atoms']):
        symbol = atom_record['element']
        element_row = elements[symbol]
        atom = add_sphere(
            f"Custom_{uuid.UUID(record['id']).hex[:8]}_Atom_{index + 1}_{symbol}",
            radii.get(symbol, 0.25),
            atom_record['position'],
            collection,
            materials[symbol],
        )
        atom['element_symbol'] = symbol
        atom['atomic_number'] = element_row[0]
        atom['element_name'] = element_row[2]
        atom['formal_charge'] = atom_record['charge']
        atom['atomlab_custom_molecule_id'] = record['id']
        atom['atomlab_custom_atom_index'] = index
        atom['atomlab_bond_order_sum'] = sum(
            bond['order'] for bond in record['bonds'] if index in (bond['a'], bond['b'])
        )
        atom['molecule_formula'] = record['formula']
        atom_objects.append(atom)

    for bond_index, bond in enumerate(record['bonds'], start=1):
        first, second, order = bond['a'], bond['b'], bond['order']
        start = Vector(record['atoms'][first]['position'])
        end = Vector(record['atoms'][second]['position'])
        vector = end - start
        length = vector.length
        if length < 0.05:
            raise ValueError(f"Saved bond {first + 1}-{second + 1} has zero/near-zero length.")
        direction = vector.normalized()
        perpendicular = _stable_perpendicular(direction)
        offsets = {1: (0.0,), 2: (-0.045, 0.045), 3: (-0.07, 0.0, 0.07)}[order]
        for strand, offset in enumerate(offsets, start=1):
            bpy.ops.mesh.primitive_cylinder_add(
                vertices=16, radius=0.032, depth=length,
                location=tuple((start + end) * 0.5 + perpendicular * offset),
            )
            rod = bpy.context.object
            rod.name = f"Custom_{uuid.UUID(record['id']).hex[:8]}_Bond_{bond_index}_{strand}"
            rod.rotation_euler = direction.to_track_quat('Z', 'Y').to_euler()
            rod.data.materials.append(bond_material)
            for polygon in rod.data.polygons:
                polygon.use_smooth = True
            put_in_collection(rod, collection)
            rod['atomlab_custom_molecule_id'] = record['id']
            rod['atom_a_name'] = atom_objects[first].name
            rod['atom_b_name'] = atom_objects[second].name
            rod['bond_order'] = order
            rod['bond_length_blender_units'] = length

    scene['atomlab_model'] = "User-saved molecule"
    scene['atomlab_custom_molecule_id'] = record['id']
    scene['atomlab_molecule_name'] = record['name']
    scene['atomlab_formula'] = record['formula']
    scene['atomlab_custom_molecule_net_charge'] = record['net_charge']
    return record


def _unit_vector(vector, fallback=Vector((1.0, 0.0, 0.0))):
    return vector.normalized() if vector.length > 1e-8 else fallback.copy()


def _stable_perpendicular(axis):
    reference = Vector((0.0, 0.0, 1.0))
    if abs(axis.dot(reference)) > 0.9:
        reference = Vector((0.0, 1.0, 0.0))
    return _unit_vector(reference - axis * reference.dot(axis), Vector((1.0, 0.0, 0.0)))


def _place_conical_neighbors(center, atoms, angle_degrees, bond_length, just_added=None):
    old_atoms = [atom for atom in atoms if atom != just_added]
    old_directions = [
        _unit_vector(atom.matrix_world.translation - center.matrix_world.translation)
        for atom in old_atoms
    ]
    all_directions = [
        _unit_vector(atom.matrix_world.translation - center.matrix_world.translation)
        for atom in atoms
    ]
    source = old_directions if old_directions else all_directions
    axis = _unit_vector(sum(source, Vector((0.0, 0.0, 0.0))), Vector((0.0, 0.0, 1.0)))
    if source:
        projected = source[0] - axis * source[0].dot(axis)
        basis = _unit_vector(projected, _stable_perpendicular(axis))
    else:
        basis = _stable_perpendicular(axis)
    side = _unit_vector(axis.cross(basis), _stable_perpendicular(axis))
    theta = math.radians(angle_degrees)
    cos_polar = math.sqrt(max(0.0, (1.0 + 2.0 * math.cos(theta)) / 3.0))
    sin_polar = math.sqrt(max(0.0, 1.0 - cos_polar * cos_polar))
    directions = [
        axis * cos_polar + sin_polar * (
            basis * math.cos(2.0 * math.pi * k / 3.0)
            + side * math.sin(2.0 * math.pi * k / 3.0)
        )
        for k in range(3)
    ]
    assignments = {}
    available = set(range(3))
    for atom, old_direction in zip(old_atoms, old_directions):
        index = max(available, key=lambda i: old_direction.dot(directions[i]))
        assignments[atom] = index
        available.remove(index)
    if just_added in atoms:
        assignments[just_added] = next(iter(available))
    else:
        remaining_atoms = [atom for atom in atoms if atom not in assignments]
        for atom, index in zip(remaining_atoms, sorted(available)):
            assignments[atom] = index
    for atom, index in assignments.items():
        atom.location = center.location + directions[index] * bond_length
    return True


def _place_tetrahedral_neighbors(center, atoms, bond_length, just_added=None):
    if len(atoms) != 4:
        return False
    old_atoms = [atom for atom in atoms if atom != just_added]
    old_directions = [
        _unit_vector(atom.matrix_world.translation - center.matrix_world.translation)
        for atom in old_atoms
    ]
    if just_added in atoms and len(old_atoms) == 3:
        axis = _unit_vector(sum(old_directions, Vector((0.0, 0.0, 0.0))), Vector((0.0, 0.0, 1.0)))
        projected = old_directions[0] - axis * old_directions[0].dot(axis)
        basis = _unit_vector(projected, _stable_perpendicular(axis))
        side = _unit_vector(axis.cross(basis), _stable_perpendicular(axis))
        cos_polar = 1.0 / 3.0
        sin_polar = math.sqrt(8.0 / 9.0)
        directions = [
            axis * cos_polar + sin_polar * (
                basis * math.cos(2.0 * math.pi * k / 3.0)
                + side * math.sin(2.0 * math.pi * k / 3.0)
            )
            for k in range(3)
        ]
        available = set(range(3))
        assignments = {}
        for atom, old_direction in zip(old_atoms, old_directions):
            index = max(available, key=lambda i: old_direction.dot(directions[i]))
            assignments[atom] = directions[index]
            available.remove(index)
        assignments[just_added] = -axis
    else:
        current = _unit_vector(
            atoms[0].matrix_world.translation - center.matrix_world.translation,
            Vector((0.0, 0.0, 1.0)),
        )
        basis = _stable_perpendicular(current)
        side = _unit_vector(current.cross(basis), _stable_perpendicular(current))
        radial = math.sqrt(8.0 / 9.0)
        directions = [current] + [
            -current / 3.0 + radial * (
                basis * math.cos(2.0 * math.pi * k / 3.0)
                + side * math.sin(2.0 * math.pi * k / 3.0)
            )
            for k in range(3)
        ]
        assignments = {atoms[0]: directions[0]}
        for atom, direction in zip(atoms[1:], directions[1:]):
            assignments[atom] = direction
    for atom, direction in assignments.items():
        atom.location = center.location + direction * bond_length
    return True


def _builder_ideal_bond_length(center_symbol, neighbor, order):
    neighbor_symbol = neighbor.get("element_symbol", "")
    key = (center_symbol, neighbor_symbol, order)
    reverse_key = (neighbor_symbol, center_symbol, order)
    known = {
        ('C', 'H', 1): 1.09, ('C', 'C', 1): 1.54, ('C', 'C', 2): 1.34,
        ('C', 'C', 3): 1.20, ('C', 'O', 1): 1.43, ('C', 'O', 2): 1.16,
        ('C', 'N', 1): 1.47, ('C', 'N', 2): 1.29, ('C', 'N', 3): 1.16,
        ('O', 'H', 1): 0.9584, ('N', 'H', 1): 1.012,
        ('O', 'O', 1): 1.48, ('O', 'O', 2): 1.21,
        ('N', 'N', 1): 1.45, ('N', 'N', 2): 1.25, ('N', 'N', 3): 1.10,
    }
    length = known.get(key, known.get(reverse_key))
    if length is None:
        r1 = COVALENT_RADII_A.get(center_symbol, 0.8)
        r2 = COVALENT_RADII_A.get(neighbor_symbol, 0.8)
        length = (r1 + r2) * {1: 1.0, 2: 0.90, 3: 0.82}.get(order, 1.0)
    return length * 1.4


def _place_linear_neighbors(center, records, just_added=None):
    atoms = [atom for atom, _order in records]
    if len(atoms) != 2:
        return False
    if just_added in atoms:
        old_atom = next(atom for atom in atoms if atom != just_added)
        axis = _unit_vector(old_atom.matrix_world.translation - center.matrix_world.translation)
        old_length = next(_builder_ideal_bond_length(center.get("element_symbol"), a, o)
                          for a, o in records if a == old_atom)
        new_length = next(_builder_ideal_bond_length(center.get("element_symbol"), a, o)
                          for a, o in records if a == just_added)
        old_atom.location = center.location + axis * old_length
        just_added.location = center.location - axis * new_length
    else:
        first, second = atoms
        v1 = first.matrix_world.translation - center.matrix_world.translation
        v2 = second.matrix_world.translation - center.matrix_world.translation
        axis = _unit_vector(v1 - v2, _unit_vector(v1))
        for atom, _order in records:
            length = _builder_ideal_bond_length(center.get("element_symbol"), atom, _order)
            direction = axis if atom == first else -axis
            atom.location = center.location + direction * length
    return True


def _place_planar_neighbors(center, records, just_added=None):
    atoms = [atom for atom, _order in records]
    if len(atoms) != 3:
        return False
    old_atoms = [atom for atom in atoms if atom != just_added]
    vectors = [
        _unit_vector(atom.matrix_world.translation - center.matrix_world.translation)
        for atom in old_atoms
    ]
    anchor = vectors[0] if vectors else _unit_vector(
        atoms[0].matrix_world.translation - center.matrix_world.translation,
        Vector((1.0, 0.0, 0.0)),
    )
    if len(vectors) > 1:
        normal = _unit_vector(vectors[0].cross(vectors[1]), _stable_perpendicular(anchor))
    else:
        normal = _stable_perpendicular(anchor)
    in_plane = _unit_vector(normal.cross(anchor), _stable_perpendicular(anchor))
    angle = 2.0 * math.pi / 3.0
    directions = [
        anchor,
        anchor * math.cos(angle) + in_plane * math.sin(angle),
        anchor * math.cos(angle) - in_plane * math.sin(angle),
    ]
    assignments = {}
    available = set(range(3))
    for atom, old_direction in zip(old_atoms, vectors):
        index = max(available, key=lambda i: old_direction.dot(directions[i]))
        assignments[atom] = index
        available.remove(index)
    if just_added in atoms:
        assignments[just_added] = next(iter(available))
    else:
        for atom, index in zip((a for a in atoms if a not in assignments), sorted(available)):
            assignments[atom] = index
    for atom, index in assignments.items():
        order = next(order for candidate, order in records if candidate == atom)
        length = _builder_ideal_bond_length(center.get("element_symbol"), atom, order)
        atom.location = center.location + directions[index] * length
    return True


def optimize_common_geometry(scene, center, just_added=None):
    """Apply idealized common geometries after bond edits; not a quantum optimization."""
    neighbors_by_name = {}
    for rod in scene.objects:
        if not rod.get("atomlab_builder_bond", False):
            continue
        a = scene.objects.get(rod.get("atom_a_name", ""))
        b = scene.objects.get(rod.get("atom_b_name", ""))
        if a is None or b is None:
            continue
        if a == center:
            neighbors_by_name[b.name] = (b, int(rod.get("bond_order", 1)))
        elif b == center:
            neighbors_by_name[a.name] = (a, int(rod.get("bond_order", 1)))

    records = list(neighbors_by_name.values())
    if not records:
        return False
    atoms = [atom for atom, _order in records]
    element = center.get("element_symbol")
    charge = int(center.get("formal_charge", 0))
    total_order = sum(order for _atom, order in records)

    # Carbon: linear CO2/other sp-like center, trigonal planar sp2 center, or CH4.
    if element == 'C' and charge == 0:
        if len(records) == 2 and total_order == 4:
            return _place_linear_neighbors(center, records, just_added)
        if len(records) == 3 and total_order == 4:
            return _place_planar_neighbors(center, records, just_added)
        if (len(records) == 4 and total_order == 4
                and all(order == 1 and atom.get("element_symbol") == 'H'
                        for atom, order in records)):
            return _place_tetrahedral_neighbors(
                center, atoms, 1.09 * 1.4, just_added
            )

    # Remaining supported VSEPR templates require four/electron-pair partners H.
    if any(order != 1 or atom.get("element_symbol") != 'H' for atom, order in records):
        return False

    # Water: neutral O with two single O-H bonds, angle 104.5 degrees.
    if element == 'O' and charge == 0 and len(atoms) == 2:
        target_length = 0.9584 * 1.4
        angle = math.radians(104.5)
        if just_added in atoms:
            old_atom = next(atom for atom in atoms if atom != just_added)
            old_direction = _unit_vector(
                old_atom.matrix_world.translation - center.matrix_world.translation
            )
            perpendicular = _stable_perpendicular(old_direction)
            old_atom.location = center.location + old_direction * target_length
            new_direction = old_direction * math.cos(angle) + perpendicular * math.sin(angle)
            just_added.location = center.location + new_direction * target_length
        else:
            v1 = _unit_vector(atoms[0].matrix_world.translation - center.matrix_world.translation)
            v2 = _unit_vector(atoms[1].matrix_world.translation - center.matrix_world.translation)
            axis = _unit_vector(v1 + v2, _stable_perpendicular(v1))
            difference = v1 - v2
            side = _unit_vector(difference - axis * difference.dot(axis), _stable_perpendicular(axis))
            half = angle * 0.5
            directions = (
                axis * math.cos(half) + side * math.sin(half),
                axis * math.cos(half) - side * math.sin(half),
            )
            for atom, direction in zip(atoms, directions):
                atom.location = center.location + direction * target_length
        return True

    # Hydronium, ammonia, ammonium.
    if element == 'O' and charge == 1 and len(atoms) == 3:
        return _place_conical_neighbors(center, atoms, 113.6, 0.98 * 1.4, just_added)
    if element == 'N' and charge == 0 and len(atoms) == 3:
        return _place_conical_neighbors(center, atoms, 107.8, 1.012 * 1.4, just_added)
    if element == 'N' and charge == 1 and len(atoms) == 4:
        return _place_tetrahedral_neighbors(center, atoms, 1.03 * 1.4, just_added)
    return False

def create_builder_bond(context, bond_order):
    selected = [obj for obj in context.selected_objects
                if obj.get("atomlab_builder_atom", False)]
    if len(selected) != 2:
        raise ValueError("Select exactly two AtomLab builder atoms (Shift-click).")
    atom_a, atom_b = selected
    if atom_a == atom_b:
        raise ValueError("Choose two different atoms.")

    collection = get_or_create_builder_collection(context.scene)
    pair = tuple(sorted((atom_a.name, atom_b.name)))
    for obj in context.scene.objects:
        if obj.get("atomlab_builder_bond", False):
            existing_pair = tuple(sorted((obj.get("atom_a_name", ""), obj.get("atom_b_name", ""))))
            if existing_pair == pair:
                raise ValueError("These atoms are already bonded; delete that bond before replacing it.")

    order = int(bond_order)
    symbol_a, symbol_b = atom_a["element_symbol"], atom_b["element_symbol"]
    valence_a = int(atom_a.get("atomlab_bond_order_sum", 0))
    valence_b = int(atom_b.get("atomlab_bond_order_sum", 0))
    warnings = []
    valence_violations = []
    for atom, symbol, valence in (
        (atom_a, symbol_a, valence_a), (atom_b, symbol_b, valence_b)
    ):
        charge = int(atom.get("formal_charge", 0))
        max_valence = simple_valence_limit(symbol, charge)
        if max_valence is not None and valence + order > max_valence:
            valence_violations.append(
                f"{symbol} at charge {charge:+d} would exceed the simple valence limit {max_valence}"
            )
    if valence_violations:
        if not context.scene.atomlab_allow_valence_override:
            raise ValueError(
                "Bond blocked by the simple valence check: "
                + "; ".join(valence_violations)
                + ". Set the correct formal charge or enable the explicit override."
            )
        warnings.extend(valence_violations)

    start = atom_a.matrix_world.translation.copy()
    end = atom_b.matrix_world.translation.copy()
    vector = end - start
    length = vector.length
    if length < 0.05:
        raise ValueError("Atoms are too close or overlap; move them apart before bonding.")
    direction = vector.normalized()

    radii_a = COVALENT_RADII_A.get(symbol_a)
    radii_b = COVALENT_RADII_A.get(symbol_b)
    if radii_a is not None and radii_b is not None:
        expected = (radii_a + radii_b) * 1.4
        if length < expected * 0.4 or length > expected * 2.5:
            warnings.append("bond length is far from the rough covalent-radius estimate")

    reference = Vector((0.0, 0.0, 1.0))
    if abs(direction.dot(reference)) > 0.90:
        reference = Vector((0.0, 1.0, 0.0))
    perpendicular = direction.cross(reference).normalized()
    offsets = {1: (0.0,), 2: (-0.045, 0.045), 3: (-0.07, 0.0, 0.07)}[order]
    material = get_or_create_emission_material(
        "AtomLab_Builder_Bond_Mat", (0.62, 0.72, 0.86), 0.55
    )
    for strand, offset in enumerate(offsets, start=1):
        bpy.ops.mesh.primitive_cylinder_add(
            vertices=16,
            radius=0.032,
            depth=length,
            location=tuple((start + end) * 0.5 + perpendicular * offset),
        )
        rod = bpy.context.object
        rod.name = f"Builder_Bond_{atom_a.name}_{atom_b.name}_{strand}"
        rod.rotation_euler = direction.to_track_quat('Z', 'Y').to_euler()
        rod.data.materials.append(material)
        for polygon in rod.data.polygons:
            polygon.use_smooth = True
        put_in_collection(rod, collection)
        rod["atomlab_builder_bond"] = True
        rod["atom_a_name"] = atom_a.name
        rod["atom_b_name"] = atom_b.name
        rod["bond_order"] = order
        rod["bond_length_blender_units"] = length
        rod["bond_base_length"] = length
        rod["bond_offset"] = offset

    atom_a["atomlab_bond_order_sum"] = valence_a + order
    atom_b["atomlab_bond_order_sum"] = valence_b + order
    if valence_a == 0 and valence_b > 0:
        just_added = atom_a
    elif valence_b == 0 and valence_a > 0:
        just_added = atom_b
    else:
        just_added = None
    for center in (atom_a, atom_b):
        if center.get("element_symbol") in {'O', 'N', 'C'}:
            optimize_common_geometry(context.scene, center, just_added)
    return length, warnings


_BOND_HANDLER_RUNNING = False


def atomlab_bond_update_handler(scene, depsgraph):
    """Keep preset and builder bond rods aligned with their two endpoint atoms."""
    global _BOND_HANDLER_RUNNING
    if _BOND_HANDLER_RUNNING:
        return
    _BOND_HANDLER_RUNNING = True
    try:
        for rod in scene.objects:
            if not rod.get("atomlab_builder_bond", False):
                continue
            atom_a = scene.objects.get(rod.get("atom_a_name", ""))
            atom_b = scene.objects.get(rod.get("atom_b_name", ""))
            if atom_a is None or atom_b is None:
                continue
            start = atom_a.matrix_world.translation.copy()
            end = atom_b.matrix_world.translation.copy()
            vector = end - start
            length = vector.length
            if length < 0.05:
                continue
            direction = vector.normalized()
            reference = Vector((0.0, 0.0, 1.0))
            if abs(direction.dot(reference)) > 0.90:
                reference = Vector((0.0, 1.0, 0.0))
            perpendicular = direction.cross(reference).normalized()
            offset = float(rod.get("bond_offset", 0.0))
            new_location = (start + end) * 0.5 + perpendicular * offset
            if (rod.location - new_location).length > 1e-6:
                rod.location = new_location
            new_rotation = direction.to_track_quat('Z', 'Y').to_euler()
            angle_delta = sum(
                abs((a - b + math.pi) % (2.0 * math.pi) - math.pi)
                for a, b in zip(rod.rotation_euler, new_rotation)
            )
            if angle_delta > 1e-6:
                rod.rotation_euler = new_rotation
            base_length = float(rod.get("bond_base_length", length))
            new_scale_z = length / base_length if base_length > 1e-8 else 1.0
            if abs(rod.scale.z - new_scale_z) > 1e-6:
                rod.scale.z = new_scale_z
    finally:
        _BOND_HANDLER_RUNNING = False


def remove_atomlab_bond_handler():
    for handler in tuple(bpy.app.handlers.depsgraph_update_post):
        if getattr(handler, "__name__", "") == "atomlab_bond_update_handler":
            bpy.app.handlers.depsgraph_update_post.remove(handler)


def clear_builder_collection(scene):
    collection = bpy.data.collections.get("AtomLab_Builder")
    if collection:
        for obj in tuple(collection.objects):
            bpy.data.objects.remove(obj, do_unlink=True)
        bpy.data.collections.remove(collection, do_unlink=True)


def create_nucleus_visual(collection, z, nucleus_mode, proton_mat, neutron_mat):
    atomic_number, symbol, name, mass_number = get_element(z)
    neutron_count = mass_number - atomic_number
    nucleus_extent = 0.13 + 0.023 * (mass_number ** (1.0 / 3.0))
    isotope = f"{symbol}-{mass_number}"

    if nucleus_mode == 'MARKERS':
        nucleon_radius = max(0.035, 0.065 * (12.0 / mass_number) ** (1.0 / 3.0))
        positions = nucleon_positions(mass_number, nucleus_extent)
        for index, location in enumerate(positions):
            choose_proton = ((index + 1) * atomic_number // mass_number) > (index * atomic_number // mass_number)
            if choose_proton:
                particle_name, material = "Proton", proton_mat
            else:
                particle_name, material = "Neutron", neutron_mat
            particle = add_sphere(
                f"{symbol}_{particle_name}_{index + 1}",
                nucleon_radius,
                location,
                collection,
                material,
            )
            particle["particle_type"] = particle_name.lower()
            particle["isotope"] = isotope
            particle["educational_note"] = (
                "Illustrative placement only; not a measured nuclear configuration."
            )
    else:
        add_density_cloud(
            f"{symbol}_Proton_Density_Approx",
            atomic_number,
            nucleus_extent,
            collection,
            proton_mat,
            seed=atomic_number * 1000 + mass_number,
            isotope=isotope,
        )
        add_density_cloud(
            f"{symbol}_Neutron_Density_Approx",
            neutron_count,
            nucleus_extent,
            collection,
            neutron_mat,
            seed=atomic_number * 2000 + mass_number,
            isotope=isotope,
        )

    return atomic_number, symbol, name, mass_number, neutron_count, nucleus_extent


def create_bohr_schematic(context, z, nucleus_mode, charge=0):
    scene = context.scene
    electron_count = z - charge
    shell_counts = get_shell_occupancies(z, charge)
    collection = create_model_collection(scene)
    atomic_number, symbol, name, mass_number = get_element(z)
    neutron_count = mass_number - atomic_number
    scene.frame_start = 1
    if scene.frame_end < 250:
        scene.frame_end = 250

    # High-contrast, color-blind-conscious-by-role palette; diffuse_color also
    # makes the colors visible in Blender's Solid viewport shading.
    orbit_mat = get_or_create_emission_material(
        "AtomLab_Orbit_Mat", (0.10, 0.38, 0.82), 0.8
    )
    electron_mat = get_or_create_emission_material(
        "AtomLab_Electron_Mat", (0.05, 0.78, 1.0), 1.8
    )
    proton_mat = get_or_create_emission_material(
        "AtomLab_Proton_Mat", (1.0, 0.22, 0.12), 1.3
    )
    neutron_mat = get_or_create_emission_material(
        "AtomLab_Neutron_Mat", (0.48, 0.30, 0.95), 1.3
    )

    # The default isotope is explicit, so A = Z + N and neutron count is defined.
    # Both nucleus views are schematic; the whole model is deliberately not to scale.
    atomic_number, symbol, name, mass_number, neutron_count, nucleus_extent = (
        create_nucleus_visual(collection, z, nucleus_mode, proton_mat, neutron_mat)
    )

    shell_radii = {1: 0.85, 2: 1.35, 3: 1.9, 4: 2.5, 5: 3.15, 6: 3.8, 7: 4.45, 8: 5.1}

    for shell_index, electron_count in shell_counts:
        radius = shell_radii[shell_index]
        add_orbit(f"{symbol}_Shell_{shell_index}_Schematic", radius, collection, orbit_mat)

        # Qualitative speed cue, not a calculated electron speed for this atom.
        # Larger shells rotate more slowly; the arbitrary radii are not to scale.
        period_frames = 249.0 * (shell_index ** 2) / 4.0

        for e_index in range(electron_count):
            angle = (2.0 * math.pi * e_index) / electron_count
            loc = (radius * math.cos(angle), radius * math.sin(angle), 0.0)
            electron = add_sphere(
                f"{symbol}_Electron_n{shell_index}_{e_index + 1}",
                0.052,
                loc,
                collection,
                electron_mat,
            )
            animate_electron_schematic(electron, radius, angle, period_frames)

    scene["atomlab_model"] = "Bohr schematic (ion/atom)"
    scene["atomlab_element_info"] = f"{name} ({symbol}), Z={atomic_number}"
    scene["atomlab_isotope"] = f"{symbol}-{mass_number}"
    scene["atomlab_protons"] = atomic_number
    scene["atomlab_neutrons"] = neutron_count
    scene["atomlab_model_ionic_charge"] = charge
    scene["atomlab_electron_count"] = electron_count
    scene["atomlab_electron_configuration_approx"] = format_electron_configuration(z, charge)
    scene["atomlab_configuration_evidence"] = electron_configuration_evidence(atomic_number, charge)
    scene["atomlab_builder_net_charge"] = builder_net_charge(scene)
    return shell_counts


def create_quantum_orbital(context, z, nucleus_mode, orbital_id):
    scene = context.scene
    charge = scene.atomlab_ion_charge
    subshells = get_subshell_occupancies(z, charge)
    n, l, orientation, orbital_name = ORBITAL_DEFINITIONS[orbital_id]
    occupancy = subshells.get((n, l), 0)
    allow_unoccupied = scene.atomlab_show_unoccupied
    if occupancy == 0 and not allow_unoccupied:
        raise ValueError(
            f"{orbital_name} is empty for this atom/ion configuration; "
            "enable 'Show unoccupied orbitals' to visualize its shape."
        )
    collection = create_model_collection(scene)
    atomic_number, symbol, name, mass_number = get_element(z)
    neutron_count = mass_number - atomic_number

    proton_mat = get_or_create_emission_material(
        "AtomLab_Proton_Mat", (1.0, 0.22, 0.12), 1.3
    )
    neutron_mat = get_or_create_emission_material(
        "AtomLab_Neutron_Mat", (0.48, 0.30, 0.95), 1.3
    )
    electron_mat = get_or_create_emission_material(
        "AtomLab_Electron_Mat", (0.05, 0.78, 1.0), 1.6
    )
    create_nucleus_visual(collection, z, nucleus_mode, proton_mat, neutron_mat)

    zeff = slater_effective_charge(
        atomic_number, subshells, n, l, target_electron_included=(occupancy > 0)
    )
    points = sample_hydrogenic_orbital(
        n=n,
        l=l,
        orientation=orientation,
        zeff=zeff,
        count=1400,
        seed=atomic_number * 100 + n * 10 + l + (sum(ord(ch) for ch in orientation) if orientation else 0),
    )
    cloud = add_orbital_point_mesh(
        f"{symbol}_{orbital_name}_Probability_Density_Approx",
        points,
        collection,
        electron_mat,
        orbital_name,
        occupancy,
        zeff,
    )
    cloud["element"] = f"{name} ({symbol})"
    cloud["isotope"] = f"{symbol}-{mass_number}"
    cloud["unoccupied_in_ground_state"] = (occupancy == 0)
    cloud["ionic_charge"] = charge
    cloud["electron_count"] = atomic_number - charge
    cloud["subshell"] = f"{n}{'spdf'[l]}"
    cloud["subshell_capacity"] = 2 * (2 * l + 1)
    cloud["electron_configuration_approx"] = format_electron_configuration(atomic_number, charge)
    cloud["configuration_evidence_note"] = electron_configuration_evidence(atomic_number, charge)
    if occupancy == 0:
        cloud["educational_note"] = (
            "Unoccupied orbital shape for teaching; not the ground-state electron density."
        )

    scene["atomlab_model"] = "Approximate quantum orbital probability density"
    scene["atomlab_element_info"] = f"{name} ({symbol}), Z={atomic_number}"
    scene["atomlab_isotope"] = f"{symbol}-{mass_number}"
    scene["atomlab_protons"] = atomic_number
    scene["atomlab_neutrons"] = neutron_count
    scene["atomlab_orbital_info"] = orbital_name
    scene["atomlab_model_ionic_charge"] = charge
    scene["atomlab_electron_count"] = atomic_number - charge
    scene["atomlab_electron_configuration_approx"] = format_electron_configuration(atomic_number, charge)
    scene["atomlab_configuration_evidence"] = electron_configuration_evidence(atomic_number, charge)
    scene["atomlab_builder_net_charge"] = builder_net_charge(scene)
    return orbital_name, occupancy, zeff


class ATOMLAB_OT_create_bohr(bpy.types.Operator):
    bl_idname = "atomlab.create_bohr"
    bl_label = "Create Bohr Schematic"
    bl_description = "Create a schematic animated shell model for a neutral atom or approximate ion"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        z = int(context.scene.atomlab_element)
        nucleus_mode = context.scene.atomlab_nucleus_mode
        charge = context.scene.atomlab_ion_charge
        try:
            shell_counts = create_bohr_schematic(context, z, nucleus_mode, charge)
        except Exception as exc:
            self.report({'ERROR'}, f"AtomLab could not create the model: {exc}")
            return {'CANCELLED'}

        atomic_number, symbol, name, mass_number = get_element(z)
        shell_text = ", ".join(f"n={n}: {count}" for n, count in shell_counts)
        charge_label = f"^{charge:+d}" if charge else ""
        self.report(
            {'INFO'},
            f"Created {symbol}-{mass_number}{charge_label}: {atomic_number} protons, "
            f"{mass_number - atomic_number} neutrons, {atomic_number - charge} electrons. "
            f"Shells: {shell_text}. Ion configuration is approximate; not to scale.",
        )
        # In Solid mode, use each material's diffuse color so the palette is visible.
        if context.area and context.area.type == 'VIEW_3D':
            try:
                context.space_data.shading.color_type = 'MATERIAL'
            except (AttributeError, TypeError):
                pass
        return {'FINISHED'}


class ATOMLAB_OT_create_quantum(bpy.types.Operator):
    bl_idname = "atomlab.create_quantum_orbital"
    bl_label = "Create Quantum Orbital Density"
    bl_description = "Show a sampled probability-density approximation for one occupied orbital"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        z = int(context.scene.atomlab_element)
        nucleus_mode = context.scene.atomlab_nucleus_mode
        orbital_id = context.scene.atomlab_orbital
        charge = context.scene.atomlab_ion_charge
        try:
            orbital_name, occupancy, zeff = create_quantum_orbital(
                context, z, nucleus_mode, orbital_id
            )
        except Exception as exc:
            self.report({'ERROR'}, f"Quantum orbital not created: {exc}")
            return {'CANCELLED'}

        _, symbol, _, _ = get_element(z)
        charge_label = f"^{charge:+d}" if charge else ""
        if occupancy == 0:
            message = (
                f"Displayed {orbital_name} as an unoccupied orbital shape; "
                "it is not part of this atom/ion's ground-state electron density."
            )
        else:
            message = (
                f"Created approximate {orbital_name} density for {symbol}{charge_label} "
                f"({z - charge} electrons); the subshell has {occupancy} electron(s), "
                f"Zeff≈{zeff:.2f}. "
                "No electron path is shown."
            )
        self.report({'INFO'}, message)
        if context.area and context.area.type == 'VIEW_3D':
            try:
                context.space_data.shading.color_type = 'MATERIAL'
            except (AttributeError, TypeError):
                pass
        return {'FINISHED'}


class ATOMLAB_OT_create_ionic_structure(bpy.types.Operator):
    bl_idname = "atomlab.create_ionic_structure"
    bl_label = "Create Ionic Crystal Fragment"
    bl_description = "Create the selected finite ionic crystal fragment; not a discrete molecule"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        structure_id = context.scene.atomlab_ionic_structure
        try:
            result = create_ionic_structure(context, structure_id)
        except Exception as exc:
            self.report({'ERROR'}, f"Ionic structure not created: {exc}")
            return {'CANCELLED'}
        counts_text = ", ".join(
            f"{element}: {count}" for element, count in sorted(result["ion_counts"].items())
        )
        self.report(
            {'INFO'},
            f"Created {result['formula']} {result['structure']} fragment "
            f"({counts_text}); extended lattice, not discrete molecules.",
        )
        if context.area and context.area.type == 'VIEW_3D':
            try:
                context.space_data.shading.color_type = 'MATERIAL'
            except (AttributeError, TypeError):
                pass
        return {'FINISHED'}


class ATOMLAB_OT_create_molecule(bpy.types.Operator):
    bl_idname = "atomlab.create_molecule"
    bl_label = "Create Molecule (Ball-and-Stick)"
    bl_description = "Create a curated molecule with approximate experimental geometry"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        molecule_id = context.scene.atomlab_molecule
        try:
            data = create_molecule_model(context, molecule_id)
        except Exception as exc:
            self.report({'ERROR'}, f"Molecule not created: {exc}")
            return {'CANCELLED'}
        self.report(
            {'INFO'},
            f"Created {data['formula']}: {data['geometry']}; "
            f"bond length ≈ {data['bond_length_A']} Å. Ball-and-stick is schematic.",
        )
        if context.area and context.area.type == 'VIEW_3D':
            try:
                context.space_data.shading.color_type = 'MATERIAL'
            except (AttributeError, TypeError):
                pass
        return {'FINISHED'}


class ATOMLAB_OT_prepare_builder(bpy.types.Operator):
    bl_idname = "atomlab.prepare_builder"
    bl_label = "Start / Reset Builder"
    bl_description = "Continue the current molecule or start a clean builder scene"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        mode = context.scene.atomlab_builder_start_mode
        try:
            prepare_builder(context, mode)
        except Exception as exc:
            self.report({'ERROR'}, f"Could not prepare builder: {exc}")
            return {'CANCELLED'}
        if mode == 'CLEAR':
            self.report({'INFO'}, "Started fresh; previous AtomLab model and builder objects cleared.")
        else:
            self.report({'INFO'}, "Builder ready; current preset/model was preserved.")
        return {'FINISHED'}


class ATOMLAB_OT_add_builder_atom(bpy.types.Operator):
    bl_idname = "atomlab.add_builder_atom"
    bl_label = "Add Builder Atom"
    bl_description = "Add the selected element to the interactive molecule-building area"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        try:
            atom = add_builder_atom(
                context,
                int(context.scene.atomlab_builder_element),
                int(context.scene.atomlab_builder_formal_charge),
            )
        except Exception as exc:
            self.report({'ERROR'}, f"Could not add atom: {exc}")
            return {'CANCELLED'}
        self.report(
            {'INFO'},
            f"Added {atom['element_symbol']} with formal charge "
            f"{int(atom.get('formal_charge', 0)):+d}; move it with G before bonding.",
        )
        # Reset the entry field to neutral to avoid carrying an ion charge
        # accidentally onto the next atom added.
        context.scene.atomlab_builder_formal_charge = 0
        return {'FINISHED'}


class ATOMLAB_OT_set_atom_charge(bpy.types.Operator):
    bl_idname = "atomlab.set_atom_charge"
    bl_label = "Set Charge on Selected Atom"
    bl_description = "Apply the builder formal charge value to exactly one selected AtomLab atom"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        try:
            atom, total_charge, warnings = set_selected_atom_charge(
                context, context.scene.atomlab_builder_formal_charge
            )
        except Exception as exc:
            self.report({'ERROR'}, f"Could not set formal charge: {exc}")
            return {'CANCELLED'}
        charge = int(atom.get("formal_charge", 0))
        message = f"{atom.name} charge set to {charge:+d}; builder net charge {total_charge:+d}."
        if warnings:
            self.report({'WARNING'}, message + " " + "; ".join(warnings) + ".")
        else:
            self.report({'INFO'}, message)
        return {'FINISHED'}


class ATOMLAB_OT_bond_selected_atoms(bpy.types.Operator):
    bl_idname = "atomlab.bond_selected_atoms"
    bl_label = "Bond Selected Atoms"
    bl_description = "Create a schematic bond between exactly two selected builder atoms"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        try:
            length, warnings = create_builder_bond(
                context, context.scene.atomlab_builder_bond_order
            )
        except Exception as exc:
            self.report({'ERROR'}, f"Could not create bond: {exc}")
            return {'CANCELLED'}
        message = f"Created bond; displayed length {length:.2f} Blender units."
        if warnings:
            self.report({'WARNING'}, message + " " + "; ".join(warnings) + ".")
        else:
            self.report({'INFO'}, message)
        return {'FINISHED'}


class ATOMLAB_OT_optimize_selected_geometry(bpy.types.Operator):
    bl_idname = "atomlab.optimize_selected_geometry"
    bl_label = "Optimize Selected Geometry"
    bl_description = "Apply idealized geometry for supported water, hydronium, ammonia or tetrahedral centers"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        selected = [obj for obj in context.selected_objects
                    if obj.get("atomlab_builder_atom", False)]
        supported = {'O', 'N', 'C'}
        if len(selected) != 1 or selected[0].get("element_symbol") not in supported:
            self.report({'ERROR'}, "Select one supported O, N or C center atom.")
            return {'CANCELLED'}
        if not optimize_common_geometry(context.scene, selected[0]):
            self.report(
                {'ERROR'},
                "No supported H2O, H3O+, NH3, NH4+, CH4, CO2 or sp2-carbon geometry found around this atom.",
            )
            return {'CANCELLED'}
        self.report({'INFO'}, "Applied idealized common-molecule geometry and bond lengths.")
        return {'FINISHED'}


class ATOMLAB_OT_clear_builder(bpy.types.Operator):
    bl_idname = "atomlab.clear_builder"
    bl_label = "Clear Builder"
    bl_description = "Remove only objects in the AtomLab Builder collection"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        clear_builder_collection(context.scene)
        context.scene["atomlab_builder_net_charge"] = builder_net_charge(context.scene)
        self.report({'INFO'}, "Cleared the AtomLab Builder collection.")
        return {'FINISHED'}


class ATOMLAB_OT_save_custom_molecule(bpy.types.Operator):
    bl_idname = "atomlab.save_custom_molecule"
    bl_label = "Save Builder Molecule to Library"
    bl_description = "Save molecule atoms, charges, coordinates and bond orders to both the .blend project and shared JSON library"

    def execute(self, context):
        scene = context.scene
        try:
            record = custom_record_from_builder(
                context,
                scene.atomlab_custom_molecule_name,
                scene.atomlab_custom_molecule_note,
            )
            path = save_custom_library_record_to_both(
                context, record, scene.atomlab_custom_library_replace
            )
            enum_key = 'CUSTOM_' + uuid.UUID(record['id']).hex.upper()
            if enum_key in CUSTOM_MOLECULES_BY_ID:
                scene.atomlab_saved_custom_molecule = enum_key
        except Exception as exc:
            self.report({'ERROR'}, f"Custom molecule was not saved: {exc}")
            return {'CANCELLED'}
        self.report(
            {'INFO'},
            f"Saved {record['name']} ({record['formula']}, net {record['net_charge']:+d}) "
            "to the shared JSON library and this .blend file. Save the .blend to keep its project copy.",
        )
        return {'FINISHED'}


class ATOMLAB_OT_create_saved_custom_molecule(bpy.types.Operator):
    bl_idname = "atomlab.create_saved_custom_molecule"
    bl_label = "Create Saved Molecule"
    bl_description = "Build the selected molecule from the custom library"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        enum_key = context.scene.atomlab_saved_custom_molecule
        record = CUSTOM_MOLECULES_BY_ID.get(enum_key)
        if record is None:
            try:
                refresh_custom_molecule_catalog()
                record = CUSTOM_MOLECULES_BY_ID.get(enum_key)
            except Exception as exc:
                self.report({'ERROR'}, f"Could not read custom-molecule libraries: {exc}")
                return {'CANCELLED'}
        if record is None:
            self.report({'ERROR'}, "Select a saved molecule first.")
            return {'CANCELLED'}
        try:
            create_saved_custom_molecule(context, record)
        except Exception as exc:
            self.report({'ERROR'}, f"Saved molecule was not created: {exc}")
            return {'CANCELLED'}
        self.report({'INFO'}, f"Created {record['name']} ({record['formula']}); net charge {record['net_charge']:+d}.")
        return {'FINISHED'}


class ATOMLAB_OT_import_custom_molecule_library(bpy.types.Operator):
    bl_idname = "atomlab.import_custom_molecule_library"
    bl_label = "Import Molecule Library JSON"
    bl_description = "Merge a validated AtomLab JSON library into both project and shared storage"
    bl_options = {'REGISTER'}

    filepath: StringProperty(name="JSON Library", subtype='FILE_PATH')
    filter_glob: StringProperty(default="*.json", options={'HIDDEN'})
    replace_conflicts: BoolProperty(
        name="Replace name/id conflicts",
        description="Explicitly replace custom records with the same id or name",
        default=False,
    )

    def invoke(self, context, _event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

    def execute(self, context):
        try:
            with open(bpy.path.abspath(self.filepath), 'r', encoding='utf-8') as handle:
                payload = json.load(handle)
            path, count = import_custom_library_to_both(payload, self.replace_conflicts)
        except Exception as exc:
            self.report({'ERROR'}, f"JSON library was not imported: {exc}")
            return {'CANCELLED'}
        self.report({'INFO'}, f"Imported/merged {count} molecule record(s); shared library: {path}. Save the .blend to persist project copies.")
        return {'FINISHED'}


class ATOMLAB_OT_export_custom_molecule_library(bpy.types.Operator):
    bl_idname = "atomlab.export_custom_molecule_library"
    bl_label = "Export Molecule Library JSON"
    bl_description = "Export the merged project and shared custom-molecule libraries"
    bl_options = {'REGISTER'}

    filepath: StringProperty(name="JSON Library", subtype='FILE_PATH')
    filter_glob: StringProperty(default="*.json", options={'HIDDEN'})
    filename_ext = ".json"

    def invoke(self, context, _event):
        self.filepath = os.path.join(os.path.dirname(custom_library_user_path()), "AtomLab_Custom_Molecules_Export.json")
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

    def execute(self, context):
        try:
            payload = custom_library_payload_for_export()
            if not payload['molecules']:
                raise ValueError("There are no custom molecules to export.")
            path = bpy.path.abspath(self.filepath)
            with open(path, 'w', encoding='utf-8') as handle:
                json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
                handle.write('\n')
        except Exception as exc:
            self.report({'ERROR'}, f"JSON library was not exported: {exc}")
            return {'CANCELLED'}
        self.report({'INFO'}, f"Exported {len(payload['molecules'])} molecule(s) to {path}.")
        return {'FINISHED'}


class ATOMLAB_PT_main(bpy.types.Panel):
    bl_label = "AtomLab"
    bl_idname = "ATOMLAB_PT_main"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "AtomLab"

    def draw(self, context):
        layout = self.layout
        scene = context.scene
        layout.label(text="Educational 3D atom models")
        layout.prop(scene, "atomlab_element", text="Element")
        layout.prop(scene, "atomlab_ion_charge", text="Ion charge")
        layout.prop(scene, "atomlab_nucleus_mode", text="Nucleus view")
        layout.prop(scene, "atomlab_orbital", text="Quantum orbital")
        layout.prop(scene, "atomlab_show_unoccupied", text="Show unoccupied orbitals")
        z = int(scene.atomlab_element)
        _, symbol, name, mass_number = get_element(z)
        charge = scene.atomlab_ion_charge
        try:
            ionic_subshells = get_subshell_occupancies(z, charge)
            ion_configuration = format_electron_configuration(z, charge)
            charge_error = None
        except ValueError as exc:
            ionic_subshells = get_subshell_occupancies(z)
            ion_configuration = format_electron_configuration(z)
            charge_error = str(exc)
        charge_label = f"{charge:+d}" if charge else "neutral"
        layout.label(text=f"Default isotope: {symbol}-{mass_number}")
        layout.label(
            text=f"Protons {z} | Neutrons {mass_number - z} | Electrons {z - charge} | q={charge_label}",
            icon='INFO',
        )
        configuration_box = layout.box()
        configuration_box.label(text="Electron configuration (approx.)")
        if charge_error:
            configuration_box.label(text=f"Invalid ion charge: {charge_error}", icon='ERROR')
        else:
            configuration_box.label(
                text=electron_configuration_evidence(z, charge), icon='INFO'
            )
        config_line = ""
        for term in ion_configuration.split():
            candidate = f"{config_line} {term}".strip()
            if len(candidate) > 34 and config_line:
                configuration_box.label(text=config_line)
                config_line = term
            else:
                config_line = candidate
        if config_line:
            configuration_box.label(text=config_line)
        elif charge_error is None:
            configuration_box.label(text="No electrons (bare nucleus).")

        legend = layout.box()
        legend.label(text="Color legend (schematic)")
        legend.label(text="Electron: cyan")
        legend.label(text="Proton: coral red")
        legend.label(text="Neutron: violet")
        legend.label(text="Orbit guide: blue")

        box = layout.box()
        box.label(text="Bohr representation (schematic)", icon='INFO')
        bohr_row = box.row()
        bohr_row.enabled = charge_error is None
        bohr_row.operator(ATOMLAB_OT_create_bohr.bl_idname, icon='OUTLINER_OB_MESH')

        orbital_id = scene.atomlab_orbital
        n, l, orientation, orbital_name = ORBITAL_DEFINITIONS[orbital_id]
        occupied = ionic_subshells.get((n, l), 0)
        subshell_capacity = 2 * (2 * l + 1)
        layout.label(
            text=f"Selected subshell {n}{'spdf'[l]}: {occupied}/{subshell_capacity} electrons"
        )
        quantum_row = layout.row()
        quantum_row.enabled = (
            charge_error is None and (occupied > 0 or scene.atomlab_show_unoccupied)
        )
        quantum_row.operator(ATOMLAB_OT_create_quantum.bl_idname, icon='SHADERFX')
        if occupied == 0 and scene.atomlab_show_unoccupied:
            layout.label(text="Empty orbital: shape only, not ground-state density.", icon='INFO')
        elif occupied == 0:
            layout.label(text=f"{orbital_name} is empty in the ground state.", icon='INFO')

        box = layout.box()
        box.label(text="Quantum view: sampled probability density", icon='INFO')
        box.label(text="Hydrogen-like orbital; approximate Zeff/configuration.")
        box.label(text="One orbital at a time; no electron path.")
        box.label(text="d: five real shapes; f: seven real shapes.")
        box.label(text="d/f orientations are representative basis shapes.")
        box.label(text="Electron count is for the subshell, not one orientation.")
        box.separator()
        box.label(text="Molecule library: ball-and-stick", icon='OUTLINER_OB_MESH')
        box.prop(scene, "atomlab_molecule", text="Preset")
        box.operator(ATOMLAB_OT_create_molecule.bl_idname, icon='MESH_UVSPHERE')
        box.label(text="Geometry and bonds are schematic.", icon='INFO')

        ionic_box = layout.box()
        ionic_box.label(text="Ionic structure prototypes (not molecules)", icon='INFO')
        ionic_box.prop(scene, "atomlab_ionic_structure", text="Structure")
        ionic_box.operator(
            ATOMLAB_OT_create_ionic_structure.bl_idname,
            text="Create Lattice Fragment",
            icon='MESH_UVSPHERE',
        )
        ionic_data = IONIC_STRUCTURE_DATA[scene.atomlab_ionic_structure]
        ionic_box.label(text=f"{ionic_data['formula']}: extended {ionic_data['prototype']} lattice.")
        ionic_box.label(text=f"Bulk coordination: {ionic_data['coordination']}.")
        ionic_box.label(text="Finite-fragment surface ions have fewer neighbors.")
        ionic_box.label(text="Formal charges/color-coded ions; no covalent rods.", icon='INFO')
        ionic_box.label(text="Some prototypes have mixed ionic-covalent bonding; charges are formal.")

        builder = layout.box()
        builder.label(text="Interactive molecule builder", icon='MOD_BUILD')
        builder.prop(scene, "atomlab_builder_start_mode", text="Start mode")
        builder.operator(ATOMLAB_OT_prepare_builder.bl_idname, icon='FILE_REFRESH')
        builder.prop(scene, "atomlab_builder_element", text="New atom")
        builder.prop(scene, "atomlab_builder_formal_charge", text="Formal charge")
        builder.operator(ATOMLAB_OT_add_builder_atom.bl_idname, icon='ADD')
        builder.label(text="Move atoms with G before bonding.")
        builder.label(text="Shift-select one atom to set its charge.")
        builder.operator(ATOMLAB_OT_set_atom_charge.bl_idname, icon='MODIFIER')
        builder.label(text="Examples: O +1 in H₃O⁺; N +1 in NH₄⁺.")
        builder.label(text="Over-limit bonds are blocked unless overridden.", icon='INFO')
        builder.label(text=f"Builder net charge: {int(scene.get('atomlab_builder_net_charge', 0)):+d}")
        builder.label(text="Shift-select exactly two atoms to bond.")
        builder.prop(scene, "atomlab_builder_bond_order", text="Bond order")
        builder.prop(scene, "atomlab_allow_valence_override", text="Allow unusual valence")
        builder.operator(ATOMLAB_OT_bond_selected_atoms.bl_idname, icon='CONSTRAINT')
        builder.operator(ATOMLAB_OT_optimize_selected_geometry.bl_idname, icon='MOD_SMOOTH')
        builder.operator(ATOMLAB_OT_clear_builder.bl_idname, icon='TRASH')
        builder.label(text="Bond rods follow moved atoms.", icon='INFO')
        builder.label(text="Common water, ammonia, methane and CO2 geometries optimize.", icon='INFO')

        custom_library = layout.box()
        custom_library.label(text="My Molecule Library", icon='FILE_FOLDER')
        custom_library.prop(scene, "atomlab_custom_molecule_name", text="Name")
        custom_library.label(
            text="Name preview: " + format_library_display_text(scene.atomlab_custom_molecule_name or "—")
        )
        custom_library.prop(scene, "atomlab_custom_molecule_note", text="Note")
        if scene.atomlab_custom_molecule_note:
            custom_library.label(
                text="Note preview: " + format_library_display_text(scene.atomlab_custom_molecule_note)[:120]
            )
        custom_library.label(text="Markup: _2 subscript, ^2- superscript (braces also work).", icon='INFO')
        custom_library.prop(scene, "atomlab_custom_library_replace", text="Allow explicit replacement")
        custom_library.operator(ATOMLAB_OT_save_custom_molecule.bl_idname, icon='FILE_TICK')
        custom_library.prop(scene, "atomlab_saved_custom_molecule", text="Saved molecule")
        custom_library.operator(ATOMLAB_OT_create_saved_custom_molecule.bl_idname, icon='MESH_UVSPHERE')
        io_row = custom_library.row(align=True)
        io_row.operator(ATOMLAB_OT_import_custom_molecule_library.bl_idname, text="Import JSON", icon='IMPORT')
        io_row.operator(ATOMLAB_OT_export_custom_molecule_library.bl_idname, text="Export JSON", icon='EXPORT')
        custom_library.label(text="Saved to this .blend + shared system JSON.", icon='INFO')
        custom_library.label(text="Save the .blend file to keep its embedded copy.")

        layout.separator()
        layout.label(text="Play timeline for Bohr animation; rotate view for quantum.", icon='PLAY')
        layout.label(text="Bohr shell speeds are illustrative only.", icon='INFO')
        layout.label(text="Quantum dots sample probability, not paths.", icon='INFO')
        if scene.atomlab_nucleus_mode == 'DENSITY':
            layout.label(text="Density cloud: sample dots, not nucleon counts.", icon='INFO')
            layout.label(text="Generic approximation; not measured data.", icon='INFO')
        else:
            layout.label(text="Nucleon positions are schematic, not measured.", icon='INFO')
        layout.label(text="Sizes are not to scale.", icon='ERROR')


CLASSES = (
    ATOMLAB_OT_create_bohr,
    ATOMLAB_OT_create_quantum,
    ATOMLAB_OT_create_ionic_structure,
    ATOMLAB_OT_create_molecule,
    ATOMLAB_OT_prepare_builder,
    ATOMLAB_OT_add_builder_atom,
    ATOMLAB_OT_set_atom_charge,
    ATOMLAB_OT_bond_selected_atoms,
    ATOMLAB_OT_optimize_selected_geometry,
    ATOMLAB_OT_clear_builder,
    ATOMLAB_OT_save_custom_molecule,
    ATOMLAB_OT_create_saved_custom_molecule,
    ATOMLAB_OT_import_custom_molecule_library,
    ATOMLAB_OT_export_custom_molecule_library,
    ATOMLAB_PT_main,
)

@bpy.app.handlers.persistent
def atomlab_custom_library_load_handler(_dummy):
    try:
        refresh_custom_molecule_catalog()
    except Exception as exc:
        print(f"AtomLab custom-molecule library could not refresh: {exc}")


def remove_atomlab_custom_library_load_handler():
    handlers = bpy.app.handlers.load_post
    for handler in tuple(handlers):
        if handler == atomlab_custom_library_load_handler:
            handlers.remove(handler)


def register():
    remove_atomlab_bond_handler()
    remove_atomlab_custom_library_load_handler()
    for cls in CLASSES:
        bpy.utils.register_class(cls)
    bpy.types.Scene.atomlab_element = EnumProperty(
        name="Element",
        description="Select an element; configurations are approximate beyond common cases",
        items=ELEMENT_ITEMS,
        default='6',
    )
    bpy.types.Scene.atomlab_ion_charge = IntProperty(
        name="Ion charge",
        description="Charge in elementary-charge units; ionic configurations use simplified rules",
        default=0,
        min=-3,
        max=8,
    )
    bpy.types.Scene.atomlab_nucleus_mode = EnumProperty(
        name="Nucleus display",
        description="Choose schematic nucleon markers or an approximate density cloud",
        items=NUCLEUS_MODE_ITEMS,
        default='MARKERS',
    )
    bpy.types.Scene.atomlab_orbital = EnumProperty(
        name="Quantum orbital",
        description="Select an s, p, d or f orbital; optionally include unoccupied orbitals",
        items=ORBITAL_ITEMS,
        default='1S',
    )
    bpy.types.Scene.atomlab_show_unoccupied = BoolProperty(
        name="Show unoccupied orbitals",
        description="Visualize an empty/excited-state orbital without claiming it is occupied",
        default=False,
    )
    bpy.types.Scene.atomlab_molecule = EnumProperty(
        name="Molecule preset",
        description="Choose a discrete molecule from the starter library",
        items=MOLECULE_ITEMS,
        default='H2O',
    )
    bpy.types.Scene.atomlab_ionic_structure = EnumProperty(
        name="Ionic crystal structure",
        description="Choose an extended ionic solid fragment, not a discrete molecule",
        items=IONIC_STRUCTURE_ITEMS,
        default='NACL_ROCKSALT',
    )
    bpy.types.Scene.atomlab_builder_start_mode = EnumProperty(
        name="Builder start mode",
        description="Continue the visible model or start with a clean AtomLab workspace",
        items=BUILDER_START_ITEMS,
        default='CONTINUE',
    )
    bpy.types.Scene.atomlab_builder_element = EnumProperty(
        name="Builder element",
        description="Choose an element to add to the interactive builder",
        items=ELEMENT_ITEMS,
        default='6',
    )
    bpy.types.Scene.atomlab_builder_bond_order = EnumProperty(
        name="Bond order",
        description="Choose a single, double or triple schematic bond",
        items=BOND_ORDER_ITEMS,
        default='1',
    )
    bpy.types.Scene.atomlab_builder_formal_charge = IntProperty(
        name="Formal charge",
        description="Formal charge assigned to a newly added or selected builder atom",
        default=0,
        min=-3,
        max=3,
    )
    bpy.types.Scene.atomlab_allow_valence_override = BoolProperty(
        name="Allow unusual valence",
        description="Allow bonds beyond simple valence limits; use only for advanced cases",
        default=False,
    )
    bpy.types.Scene.atomlab_custom_molecule_name = StringProperty(
        name="Custom molecule name",
        description="Name to store; use _2 for subscript and ^2- for superscript, for example H_2O or SO_4^2-",
        default="",
        maxlen=80,
    )
    bpy.types.Scene.atomlab_custom_molecule_note = StringProperty(
        name="Educational note",
        description="Optional note; _2 and ^2- markup is rendered in the preview",
        default="",
        maxlen=240,
    )
    bpy.types.Scene.atomlab_custom_library_replace = BoolProperty(
        name="Allow explicit replacement",
        description="Replace a saved custom molecule with the same name only when explicitly enabled",
        default=False,
    )
    bpy.types.Scene.atomlab_saved_custom_molecule = EnumProperty(
        name="Saved custom molecule",
        description="Choose a saved molecule from project-embedded and shared JSON libraries",
        items=custom_molecule_enum_items,
        default=0,
    )
    bpy.app.handlers.depsgraph_update_post.append(atomlab_bond_update_handler)
    bpy.app.handlers.load_post.append(atomlab_custom_library_load_handler)
    try:
        refresh_custom_molecule_catalog()
    except Exception as exc:
        print(f"AtomLab custom-molecule library not loaded: {exc}")


def unregister():
    remove_atomlab_bond_handler()
    remove_atomlab_custom_library_load_handler()
    if hasattr(bpy.types.Scene, "atomlab_element"):
        del bpy.types.Scene.atomlab_element
    if hasattr(bpy.types.Scene, "atomlab_ion_charge"):
        del bpy.types.Scene.atomlab_ion_charge
    if hasattr(bpy.types.Scene, "atomlab_nucleus_mode"):
        del bpy.types.Scene.atomlab_nucleus_mode
    if hasattr(bpy.types.Scene, "atomlab_orbital"):
        del bpy.types.Scene.atomlab_orbital
    if hasattr(bpy.types.Scene, "atomlab_show_unoccupied"):
        del bpy.types.Scene.atomlab_show_unoccupied
    if hasattr(bpy.types.Scene, "atomlab_molecule"):
        del bpy.types.Scene.atomlab_molecule
    if hasattr(bpy.types.Scene, "atomlab_ionic_structure"):
        del bpy.types.Scene.atomlab_ionic_structure
    if hasattr(bpy.types.Scene, "atomlab_builder_start_mode"):
        del bpy.types.Scene.atomlab_builder_start_mode
    if hasattr(bpy.types.Scene, "atomlab_builder_element"):
        del bpy.types.Scene.atomlab_builder_element
    if hasattr(bpy.types.Scene, "atomlab_builder_bond_order"):
        del bpy.types.Scene.atomlab_builder_bond_order
    if hasattr(bpy.types.Scene, "atomlab_builder_formal_charge"):
        del bpy.types.Scene.atomlab_builder_formal_charge
    if hasattr(bpy.types.Scene, "atomlab_allow_valence_override"):
        del bpy.types.Scene.atomlab_allow_valence_override
    if hasattr(bpy.types.Scene, "atomlab_custom_molecule_name"):
        del bpy.types.Scene.atomlab_custom_molecule_name
    if hasattr(bpy.types.Scene, "atomlab_custom_molecule_note"):
        del bpy.types.Scene.atomlab_custom_molecule_note
    if hasattr(bpy.types.Scene, "atomlab_custom_library_replace"):
        del bpy.types.Scene.atomlab_custom_library_replace
    if hasattr(bpy.types.Scene, "atomlab_saved_custom_molecule"):
        del bpy.types.Scene.atomlab_saved_custom_molecule
    for cls in reversed(CLASSES):
        try:
            bpy.utils.unregister_class(cls)
        except RuntimeError:
            pass


if __name__ == "__main__":
    # Safe re-run from Blender's Text Editor.
    try:
        unregister()
    except Exception:
        pass
    register()
