from ase.spacegroup import crystal
from ase.visualize import view
from ase.atoms import Atoms, Cell
import argparse
from os import path, mkdir
import logging
import numpy as np


ap = argparse.ArgumentParser(
    prog="Ti3AlCN",
    description="Titanium Aluminum Cyanide MAX and MXene crystal maker.",
    epilog="written by Ali AghaGhadiri")

lattice_ap = ap.add_argument_group("Lattice", "Manipulate lattice parameters")
lattice_ap.add_argument("-a", default=3.02, type=float,
                        help="Lattice parameter `a`")
lattice_ap.add_argument("-c", default=19.35, type=float,
                        help="Lattice parameter `c`")
lattice_ap.add_argument("--super-cell", default=[1, 1, 1], dest="size",
                        type=int, nargs=3, help="Repetition of unitcell")

phase_ap = ap.add_argument_group("Phase", "determine crystal phase.")
phase_ap.add_argument("--max", action="store_true", help="MAX phase")
phase_ap.add_argument("--mx", action="store_true", help="MX phase")
phase_ap.add_argument("--mono-layer", action="store_true", dest="mono",
                      help="Monolayer Structure")

output_ap = ap.add_argument_group("Output", "Determine the output format.")
output_ap.add_argument(
    "--format", choices=["xyz", "pdb"], default="xyz", help="Output format")
output_ap.add_argument("-o", "--out", help="Output directory",
                       type=str, dest="directory", default=".")
output_ap.add_argument("--show", default="hex",
                       choices=["hex", "cube"], help="View results")

manipulate_ap = ap.add_argument_group("Manipulation")
manipulate_ap.add_argument("--drill", type=float, help="Make vacancy.")

atomic_positions_mx_mono = [
    (0.0,     0.0,     0.5),  # Ti
    (0.33333, 0.66667, 0.66667),  # Ti
    (0.66667, 0.33333, 0.33333),  # Ti
    (0.66667, 0.33333, 0.58333),  # C
    (0.33333, 0.66667, 0.41667)  # N
]

atomic_positions_mx_multi = [
    (0.00000, 0.00000, 0.00000),  # Ti
    (0.00000, 0.00000, 0.50000),  # Ti
    (0.66667, 0.33333, 0.16667),  # Ti
    (0.33333, 0.66667, 0.66667),  # Ti
    (0.33333, 0.66667, 0.83333),  # Ti
    (0.66667, 0.33333, 0.33333),  # Ti
    (0.66667, 0.33333, 0.91667),  # C
    (0.33333, 0.66667, 0.41667),  # C
    (0.33333, 0.66667, 0.08333),  # N
    (0.66667, 0.33333, 0.58333),  # N
]

atomic_positions_max = [
    (0.00000, 0.00000, 0.00000),  # Ti
    (0.00000, 0.00000, 0.50000),  # Ti
    (0.66667, 0.33333, 0.16667),  # Ti
    (0.33333, 0.66667, 0.66667),  # Ti
    (0.33333, 0.66667, 0.83333),  # Ti
    (0.66667, 0.33333, 0.33333),  # Ti
    (0.66667, 0.33333, 0.75000),  # Al
    (0.33333, 0.66667, 0.25000),  # Al
    (0.66667, 0.33333, 0.91667),  # C
    (0.33333, 0.66667, 0.41667),  # C
    (0.33333, 0.66667, 0.08333),  # N
    (0.66667, 0.33333, 0.58333),  # N
]

def drill(atoms: Atoms, radius: float):
    # Calculate the center of the atoms object
    xc = atoms.positions[:, 0].max() / 2
    yc = atoms.positions[:, 1].max() / 2
    center = np.array([xc, yc])
    # Create a mask for atoms within the given radius from the center in xy-plane
    mask = np.linalg.norm(atoms.positions[:, :2] - center, axis=1) < radius
    # Remove the atoms within the given radius from the center
    del atoms[mask]


def shrink(crystal: Atoms) -> None:
    pos = crystal.get_positions()
    pos_min = np.min(pos[:, 2])
    crystal.translate([0, 0, -pos_min])
    new_c = np.max(pos[:, 2]) - pos_min
    cell = crystal.get_cell()
    cell[2, 2] = new_c
    crystal.set_cell(cell)


def make_orthogonal(crystal: Atoms):
    cell = crystal.get_cell(complete=True)
    cell[1, 0] = 0
    cell[2, 0] = 0
    cell[2, 1] = 0
    crystal.set_cell(cell)
    crystal.wrap()
    # return crystal


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    arguments = ap.parse_args()

    print(arguments)
    a = arguments.a  # Lattice parameter a
    c = arguments.c  # Lattice parameter c
    size = arguments.size
    cellpar = [a, a, c, 90, 90, 120]

    logging.info(
        f"A super cell size: {' x '.join(map(str,size))}")
    logging.info(
        "Unit cell parameters: a={}, b={}, c={}; α={}, β={}, γ={}.".format(*cellpar))

    mono_mx_formula = "Ti3CN"
    multi_max_formula = "Ti6Al2C2N2"
    multi_mx_formula = "Ti6C2N2"
    bucket = {}

    if arguments.mono and arguments.mx:
        logging.info("Making mono layer of MX structure.")
        crstl = crystal(mono_mx_formula, basis=atomic_positions_mx_mono, pbc=[1, 1, 0],
                        size=size, cellpar=cellpar, symprec=0.1, onduplicates="replace")
        shrink(crstl)
        bucket["Ti3CN.mono"] = crstl
        logging.info("MAX mono layer created.")

    if arguments.mx and not arguments.mono:
        logging.info("Making multi layer MX structure.")
        bucket["Ti3CN.multi"] = crystal(multi_mx_formula, basis=atomic_positions_mx_multi,
                                        size=size, cellpar=cellpar, symprec=0.1, onduplicates="replace")
        logging.info("MAX multi layer created.")

    if arguments.max:
        max_cellpar = [3.04, 3.04, 18.39, 90, 90, 120]
        logging.warning(
            "Cell parameter have been changed for MAX structure: a={}, b={}, c={}; α={}, β={}, γ={}.".format(*max_cellpar))
        bucket["Ti3AlCN.multi"] = crystal(multi_max_formula, basis=atomic_positions_max,
                                          size=size, cellpar=max_cellpar, symprec=0.1, onduplicates="replace")
        logging.info("MX multi layer created.")

    if arguments.drill:
        for c in bucket.values():
            radius = arguments.drill
            drill(c, radius)

    for f, c in bucket.items():
        if arguments.show == "cube":
            make_orthogonal(c)
        logging.info(f"Showing {f}")
        view(c)

        # if not path.isdir(arguments.directory):
        #     logging.warning(f"`{arguments.directory}` is not a directory!")
        #     logging.info(f"creating `{arguments.directory}`")
        #     mkdir(arguments.directory)

        # name = path.join(arguments.directory, f"{f}.{arguments.format}")
        # logging.info(f"Writing `{name}`.")
        # c.write(name)
