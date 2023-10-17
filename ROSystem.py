import argparse
from ase.io import read
import numpy as np
import subprocess
import tempfile
import os
from ase.visualize import view


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
ap.add_argument("-o", "--out", default="./ROSystem", help="Output file.")
ap.add_argument("-t", "--tolerance", default=2.0, help="Set tolerance.")
ap.add_argument("--show", action="store_true", help="View the packing result.")

membrane_ap.add_argument(
    "-m", "--membrane", required=True, help="Membrane file.")
piston_ap.add_argument("-p", "--piston", required=True, help="Piston file.")
mixture_ap.add_argument("-M", "--molecule", nargs=2, dest="molecules", action="append",
                        required=True, help="Mixture files and count.")
system_ap.add_argument("-s", "--separation", required=True, type=float,
                       help="Separation of the piston and membrane.")


def pack(packmol: str, membrane: str, piston: str, separation: float,
         molecules: list[str, int], output: str, tolerance=2.0):
    """insert mixture molecules using `packmol` executable."""
    inp = tempfile.NamedTemporaryFile(mode="w+", suffix="inp")
    out_type = os.path.basename(output).split(".")[-1]

    atoms = read(membrane)
    xh, yh, zl = np.amax(atoms.positions, axis=0)
    xl, yl, _ = np.amin(atoms.positions, axis=0)
    zh = separation

    inp.write("\n# packmol input\n")
    inp.write(f"tolerance {tolerance}\n")
    inp.write(f"output {output}\n")
    inp.write(f"filetype {out_type}\n\n")
    # membrane structure
    inp.writelines([
        f"structure {membrane}\n",
        "  number 1\n",
        f"  fixed 0. 0. 0. 0. 0. 0.\n",
        "end structure\n\n",
    ])
    # mixture structures
    for molecule, count in molecules:
        inp.writelines([
            f"structure {molecule}\n",
            f"  number {int(count)}\n",
            f"  inside box {xl} {yl} {zl} {xh} {yh} {zh}\n",
            "end structure\n\n",
        ])
    # piston structure
    inp.writelines([
        f"structure {piston}\n",
        "  number 1\n",
        f"  fixed 0. 0. {zh} 0. 0. 0.\n",
        "end structure\n",
    ])
    inp.seek(0)
    print(inp.read())
    inp.seek(0)

    subprocess.call(packmol, stdin=inp)
    inp.close()


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
    stdout = subprocess.check_output("vmd -nt -dispdev text -eofexit",
                                     stdin=inp, text=True)
    inp.close()
    print(stdout)


if __name__ == "__main__":
    arguments = ap.parse_args()
    print(arguments)

    pack(
        packmol=arguments.executable,
        membrane=arguments.membrane,
        piston=arguments.piston,
        separation=arguments.separation,
        molecules=arguments.molecules,
        output=arguments.out,
        tolerance=arguments.tolerance
    )

    lammps_data(arguments.out)

    if arguments.show:
        system_pdb = read(arguments.out)
        view(system_pdb)
