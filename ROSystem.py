import argparse
import os
import re
import subprocess
import tempfile
from math import ceil

import numpy as np
from ase.atoms import Atoms
from ase.io import read
from ase.visualize import view

from MyCrystal import drill, make_orthogonal
from solution import solvate

# Argument parser
ap = argparse.ArgumentParser(
    "ROSystem", description="Reverse Osmosis System builder.")

# Argument groups
membrane_ap = ap.add_argument_group(
    "Membrane", "Specify characteristics of membrane.")
piston_ap = ap.add_argument_group(
    "Piston", "Specify characteristics of piston.")
mixture_ap = ap.add_argument_group(
    "Mixture", "Specify characteristics of mixture.")
system_ap = ap.add_argument_group(
    "System", "Specify characteristics of system.")

# Arguments
ap.add_argument("-x", "--executable", default="packmol",
                help="Packmol executable path.")
ap.add_argument("-o", "--out", default="./ROSystem.pdb", help="Output file.")
ap.add_argument("-t", "--tolerance", default=2.0, help="Set tolerance.")
ap.add_argument("--show", action="store_true", help="View the packing result.")
ap.add_argument("-T", "--temperature", type=float,
                default=298.15, help="Solution temperature.")

membrane_ap.add_argument("-m", "--membrane", help="Membrane file.")
membrane_ap.add_argument("-d", "--drill", type=float,
                         help="Diameter of the membrane hole.")
piston_ap.add_argument("-p", "--piston", help="Piston file.")
mixture_ap.add_argument("-M", "--molecule", nargs=2, dest="molecules", action="append",
                        required=True, help="Mixture compounds and concentration.")
system_ap.add_argument("-s", "--separation", type=float,
                       help="Separation of the piston and membrane.")
system_ap.add_argument("-c", "--crystal", type=str,
                       help="Unit cell crystal file.")
system_ap.add_argument("--size", nargs=3, type=float,
                       help="Size of the chamber.")




def modify(data: str, settings: dict):
    """Makes modifications to `data` based on provided `settings`.
    `data` refers to `LAMMPS` data file."""
    # Define regular expression patterns
    regex_headers = r"(?:(?:\ ?)(?:\d+ (?:atom|bond|angle|dihedral|improper)(?:s| types)\n))+(?:(?:(?:\ ?)-?(?:\d+.\d* )+ (?:x|y|z)lo (?:x|y|z)hi)\n)+"
    regex_header_z = r"(?P<zhi>-?\d+.\d+)  zlo zhi"
    regex_masses = r"(?: Masses\n\n)(?: \d \d*.\d* # \w*\n)+"
    regex_section = r"# (?P<section>\w*) Coeffs\n#\n(?P<lines>(?:# (?:\d*)  (?:(?:[a-zA-Z]+-?){0,})\n)+)"
    regex_types = r"# (?P<id>\d+)  (?P<type>(?:\w+-?)*)"
    regex_atom = r"(?P<num>\d*) (?P<molecule_tag>\d*) (?P<atom_type>\d*) (?P<charge>\d*.\d*) (?P<xyz>.*) # (?P<type>\w*) (?P<residue>(?:\w*)?)"
    # Compile patterns
    headers = re.compile(regex_headers)
    header_z = re.compile(regex_header_z)
    masses = re.compile(regex_masses)
    section = re.compile(regex_section)
    coeff = re.compile(regex_types)
    atom = re.compile(regex_atom)

    coefficients = settings.get("coeffs", {})
    charges = settings.get("charges", {})
    modified_data: list[str] = ["# LAMMPS data modified.\n"]
    # Double up system size in z direction to create an empty chamber bellow membrane.
    _headers: str = headers.findall(data)[-1]  #
    zhi = header_z.search(_headers).group(1)
    _headers = _headers.splitlines()
    _headers[-1] = f" {-float(zhi):.6f} {zhi}  zlo zhi"
    _headers = "\n".join(_headers)

    modified_data.append(_headers)
    modified_data.append("\n\n")
    modified_data.extend(masses.findall(data))
    # Insert any coefficient in its corresponding section.
    for section_match in section.finditer(data):
        section_match_dict = section_match.groupdict()
        _section = section_match_dict["section"]
        _section_lines: list[str] = []
        # this section is commented out if no coeffs provided.
        _section_activated = False

        for line_match in coeff.finditer(section_match_dict["lines"]):
            line_match_dict = line_match.groupdict()
            _type = line_match_dict["type"]
            _section_coefficients = coefficients.get(_section)
            _type_coefficient = _section_coefficients.get(
                _type) if _section_coefficients is not None else None
            _coefficient = " ".join(
                map(str, _type_coefficient)) if _type_coefficient is not None else ""
            _section_activated = True if _coefficient else False
            _id = line_match_dict['id'] if _coefficient else f"# {line_match_dict['id']}"
            _line = f"{_id} {_coefficient} # {_type}\n"
            _section_lines.append(_line)

        _section_line = f"\n {_section} Coeffs\n\n" if _section_activated else f"\n# {_section} Coeffs\n\n"
        modified_data.append(_section_line)
        modified_data.extend(_section_lines)

    modified_data.append("\n Atoms # full\n\n")
    # Charge corrections
    for atom_match in atom.finditer(data):
        atom_dict = atom_match.groupdict()
        charge = charges.get(atom_dict["type"])
        atom_dict["charge"] = charge if charge is not None else 0
        modified_data.append(
            "{num} {molecule_tag} {atom_type} {charge:.6f} {xyz} # {type} {residue}\n".format(**atom_dict))
    # Insert anything after Bonds section without modifications.
    modified_data.extend(re.findall("\n Bonds\n(?:.*\n)+", data))

    return "".join(modified_data)


def pack(packmol: str, membrane: str | Atoms, piston: str | Atoms, separation: float,
         molecules: list[str, float], output: str, temperature: float = 298.15, tolerance=2.0):
    """
    Insert mixture molecules using `packmol` executable.

    Parameters:
    - packmol: Path to the packmol executable.
    - membrane: Path to the file containing the membrane Atoms object.
    - piston: Path to the file containing the piston Atoms object.
    - separation: Distance between the piston and membrane.
    - molecules: List of molecule types and their quantities.
    - output: Path to the output file.
    - tolerance: Tolerance for packing molecules (default is 2.0).

    Returns:
    None. Writes the packed system to the output file.
    """
    # parameter corrections
    directory = os.path.dirname(output)
    output_name = os.path.basename(output)
    solution_name = "solution_" + output_name
    solution_output = os.path.join(directory, solution_name)
    molecules = {k: float(v) for k, v in molecules}
    # find center of the membrane.
    if not isinstance(membrane, Atoms):
        membrane_atoms = read(membrane)
    else:
        membrane_atoms = membrane
    membrane_maximum = membrane_atoms.positions.max(axis=0)
    membrane_minimum = membrane_atoms.positions.min(axis=0)
    xc, yc, _ = (membrane_maximum + membrane_minimum) / 2
    zc = membrane_maximum[2]
    # calculate size of the system and generate a mixture.
    dx, dy, _ = membrane_maximum - membrane_minimum
    dz = separation
    solvate(dx, dy, dz, temperature, tolerance, output=solution_output,
            packmol=packmol, **molecules)
    # translate the bottom of the mixture on top of membrane.
    solution_atoms = read(solution_output)
    solution_maximum = solution_atoms.positions.max(axis=0)
    solution_minimum = solution_atoms.positions.min(axis=0)
    sx, sy, _ = (solution_maximum + solution_minimum) / 2
    sz = solution_minimum[2]
    solution_shift = (xc - sx, yc - sy, zc - sz + (tolerance))
    solution_atoms.translate(solution_shift)
    # translate the bottom of the piston on top of the mixture.
    if not isinstance(piston, Atoms):
        piston_atoms = read(piston)
    else:
        piston_atoms = piston
    piston_maximum = piston_atoms.positions.max(axis=0)
    piston_minimum = piston_atoms.positions.min(axis=0)
    px, py, _ = (piston_maximum + piston_minimum) / 2
    pz = piston_minimum[2]
    piston_shift = (xc - px, yc - py, zc - pz + separation + (2*tolerance))
    piston_atoms.translate(piston_shift)
    # combine the system
    system_atoms = membrane_atoms + solution_atoms + piston_atoms
    system_atoms.write(output)


def lammps_data(ROSystem: str, settings: dict = {}):
    ROSystem = ROSystem.replace('\\', '/')
    directory = os.path.dirname(ROSystem)
    file = os.path.basename(ROSystem)
    file, *_ = file.rpartition(".")
    lmpdata = f"{directory}/{file}.data"

    system = read(ROSystem)
    x0, y0, z0 = np.amin(system.positions, axis=0)
    x1, y1, z1 = np.amax(system.positions, axis=0)
    dx, dy, dz = (x1 - x0, y1 - y0, z1 - z0)
    *_, alpha, beta, gamma = system.get_cell_lengths_and_angles()

    inp = tempfile.NamedTemporaryFile(mode="w+", suffix="inp")
    inp.writelines([
        "package require topotools\n",
        f"mol new {ROSystem}\n",
        f"molinfo top set a {dx}\n",
        f"molinfo top set b {dy}\n",
        f"molinfo top set c {dz}\n",
        f"molinfo top set alpha {alpha}\n",
        f"molinfo top set beta {beta}\n",
        f"molinfo top set gamma {gamma}\n",
        "topo guessatom element name\n",
        "topo guessbonds\n",
        "topo guessangles\n",
        "topo guessdihedrals\n",
        "topo guessimpropers\n",
        f"topo writelammpsdata {lmpdata} full\n",
        "exit\n",
    ])
    inp.seek(0)
    # NOTE: The pipe does not exit automatically.
    #       it might be due to the stdin pipe which
    #       is waiting for more input!
    print("\npress ENTER\n")
    stdout = subprocess.check_output(["vmd", "-nt", "-dispdev", "text", "-eofexit"],
                                     stdin=inp, text=True)
    inp.close()
    print(stdout)

    with open(lmpdata) as lammps_data_file:
        modified_lammps_data_file = modify(lammps_data_file.read(), settings)
    
    with open(lmpdata, "w") as lammps_data_file:
        lammps_data_file.write(modified_lammps_data_file)


def mklyr(crystal_info: str, size: list[float], drill_radius: float = 0.0):
    crystal = read(crystal_info)
    a, b, c, *_ = crystal.get_cell_lengths_and_angles()
    A = ceil(size[0] / a)
    B = ceil(size[1] / b)
    layer = crystal * [A, B, 1]
    # TODO: Trim the edges to get the proper size.
    make_orthogonal(layer)
    if drill_radius:
        drill(layer, drill_radius)
    return layer


def mksystem(unit_cell: str, size: list[float], hole_diameter: float,
             packmol: str, output: str, tolerance: float,
             temperature: float = 298.15, **molecules):
    hole_radius = hole_diameter / 2
    piston = mklyr(unit_cell, size)
    membrane = mklyr(unit_cell, size, hole_radius)
    molecules = [[m, c] for m, c in molecules.items()]
    pack(packmol, membrane, piston, size[2],
         molecules, output, temperature, tolerance)


def conf_to_dict(file_path):
    import configparser
    config = configparser.ConfigParser(inline_comment_prefixes="#;")
    config.read(file_path)
    # Convert the ConfigParser object to a dictionary
    config_dict = {"coeffs": {}, "charges": {}}
    for section in config.sections():
        if section.lower() == "charges":
            section_dict = config_dict["charges"]
        else:
            coeff = config_dict["coeffs"]
            coeff.setdefault(section, {})
            section_dict = coeff[section]
        # Convert string values to float or list of floats
        for key, value in config.items(section):
            key = "-".join(map(str.capitalize, key.split("-")))
            try:
                # Try to convert to a single float
                section_dict[key] = float(value)
            except ValueError:
                # If that fails, split the string and convert each part to a float
                section_dict[key] = [float(v) for v in value.split()]
    return config_dict

if __name__ == "__main__":
    settings = conf_to_dict("config.ro.ini")
    arguments = ap.parse_args()
    print(arguments)

    if arguments.membrane and arguments.piston:

        pack(
            packmol=arguments.executable,
            membrane=arguments.membrane,
            piston=arguments.piston,
            separation=arguments.separation,
            molecules=arguments.molecules,
            output=arguments.out,
            temperature=arguments.temperature,
            tolerance=arguments.tolerance
        )


    elif arguments.crystal and arguments.size:
        molecules = {k: float(v) for k, v in arguments.molecules}
        mksystem(
            unit_cell=arguments.crystal,
            size=arguments.size,
            hole_diameter=arguments.drill,
            packmol=arguments.executable,
            output=arguments.out,
            tolerance=arguments.tolerance,
            temperature=arguments.temperature,
            **molecules
        )

    lammps_data(arguments.out, settings)

    if arguments.show:
        system_pdb = read(arguments.out)
        view(system_pdb)
