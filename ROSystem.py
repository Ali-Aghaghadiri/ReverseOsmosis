import argparse
from math import ceil
from ase.io import read
import numpy as np
import subprocess
import tempfile
import os
from ase.visualize import view
from MyCrystal import drill, make_orthogonal
from solution import solvate
from ase.atoms import Atoms


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
ap.add_argument("-T", "--temperature", type=float, default=298.15, help="Solution temperature.")

membrane_ap.add_argument("-m", "--membrane", help="Membrane file.")
membrane_ap.add_argument("-d", "--drill", type=float, help="Diameter of the membrane hole.")
piston_ap.add_argument("-p", "--piston", help="Piston file.")
mixture_ap.add_argument("-M", "--molecule", nargs=2, dest="molecules", action="append",
                        required=True, help="Mixture compounds and concentration.")
system_ap.add_argument("-s", "--separation", type=float, help="Separation of the piston and membrane.")
system_ap.add_argument("-c", "--crystal", type=str, help="Unit cell crystal file.")
system_ap.add_argument("--size", nargs=3, type=float, help="Size of the chamber.")


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


def lammps_data(ROSystem: str):
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
    pack(packmol, membrane, piston, size[2], molecules, output, temperature, tolerance)
    lammps_data(output)


if __name__ == "__main__":
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

        lammps_data(arguments.out)
    
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

    if arguments.show:
        system_pdb = read(arguments.out)
        view(system_pdb)
