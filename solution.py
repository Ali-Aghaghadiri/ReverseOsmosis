import os
import tempfile
import subprocess
import numpy as np


densities = np.array([
    [0.1, 999.85], [1, 999.9], [4, 999.97], [10, 999.7], [15, 999.1],
    [20, 998.21], [25, 997.05], [30, 995.65], [35, 994.03], [40, 992.22],
    [45, 990.21], [50, 988.04], [55, 985.69], [60, 983.2], [65, 980.55],
    [70, 977.76], [75, 974.84], [80, 971.79], [85, 968.61], [90, 965.31],
    [95, 961.89], [100, 958.35]
])

def num_water_molecules(x: float, y: float, z: float, temp: float = 298.15) -> int:
    """
    Estimate the number of water molecules in a volume at 1 atm pressure.
    x, y, z: volume dimension in angstroms.
    temp: temperature in Kelvin.
    """
    A = 1e-10  # Angstrom to meter convert factor
    N = 6.02214076e23  # Avogadro's number
    T = temp - 273.15  # Convert to C
    D = np.interp(T, densities[:, 0], densities[:, 1])  # Density of water at T (kg/m^3)
    M = 18.01528 * 1e-3  # Molar mass of water (kg/mol)
    V = (x * y * z) * A**3  # System volume (m^3)
    m = D * V  # Mass of water (g)
    moles_water = m / M  # moles of water in the system
    return int(moles_water * N)  # number of water molecules


def num_ions(concentration: float, water_count: int):
    """
    Estimate the number of ions in a hydrate system.

    Parameters:
    concentration (float): The concentration of the ion in molar.
    water_count (int): The number of water molecules in the system.

    Returns:
    int: The estimated number of ions in the system.
    """
    N = 6.02214076e23  # Avogadro's number
    moles_water = water_count / N  # Moles of water
    # Molar concentration of water is approximately 55.345 M
    volume_water_L = moles_water / 55.345  # Volume of water (L)
    moles_ions = concentration * volume_water_L  # moles of ions
    return int(moles_ions * N)  # Number of ions


def solvate(x: float, y: float, z: float, temperature: float = 298.15, tol: float = 2.0, directory: str = "./mixture/", output: str = "solvated.pdb", packmol: str = "packmol", **molecules) -> None:
    directory = os.path.abspath(directory)
    output = os.path.abspath(output)

    num_waters = num_water_molecules(x, y, z, temperature)

    inp = tempfile.NamedTemporaryFile(mode="w+", suffix="inp")
    inp.write("\n".join([
        f"tolerance {tol}",
        "filetype pdb",
        f"output {output}\n",
        f"structure {directory}/H2O.pdb",
        f"  number {num_waters}",
        f"  inside box 0. 0. 0. {x} {y} {z}\n"
        "end structure\n",
    ]))

    for ion, concentration in molecules.items():
        ions_count = num_ions(concentration, num_waters)
        inp.writelines([
            f"\nstructure {directory}/{ion}.pdb\n"
            f"  number {ions_count}\n"
            f"  inside box 0. 0. 0. {x} {y} {z}\n"
            "end structure\n"
        ])
    inp.seek(0)
    print(inp.read())
    inp.seek(0)
    subprocess.call(packmol, stdin=inp)
    inp.close


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser("Solution", description="Disolve molecules in water.")
    size_ap = ap.add_mutually_exclusive_group(required=True)
    soluble_ap = ap.add_mutually_exclusive_group(required=True)

    ap.add_argument("-d", "--soluble-dir", dest="directory",
                    default="./mixture/", help="Looking for solubles in this directory.")
    ap.add_argument("-t", "--tolerance", type=float,
                    default=2.0, help="Distance tolerance.")
    ap.add_argument("-o", "--out", default="./solvated.pdb",
                    help="Output pdb file. Default is ./solvated.pdb")
    ap.add_argument("-T", "--temperature", type=float, help="Solution temperature.")
    size_ap.add_argument("-s", "--size", type=float, nargs=3)
    soluble_ap.add_argument("-a", "--add", dest="solubles", nargs=2,
                            action="append", help="Specify type and concentration.")

    arguments = ap.parse_args()

    x, y, z = arguments.size
    solubles = {k: float(v) for k, v in arguments.solubles}
    print(solubles)

    packmol_input = solvate(
        x, y, z,
        arguments.temperature,
        arguments.tolerance,
        arguments.directory,
        arguments.out,
        **solubles
    )
