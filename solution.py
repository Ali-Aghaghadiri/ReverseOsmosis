import os
import tempfile


def num_water_molecules(x: float, y: float, z: float) -> int:
    """
    Estimate the number of water molecules in a volume.
    x, y, z: volume dimension in angstroms.
    """
    A = 1e-10  # Angstrom convert factor
    N = 6.02214076e23  # Avogadro's number
    D = 1000  # Density of water at room temperature (kg/m^3)
    M = 18.01528  # Molar mass of water (g/mol)
    V = (x * y * z) * A**3  # System volume (m^3)
    m = (D * V) * 1e3 # Mass of water (g)
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

def solvate(x: float, y: float, z: float, tol: float = 2.0, directory: str = "./mixture/", **kwargs):
    """Returns packmol input temp file"""
    directory = os.path.abspath(directory)
    num_waters = num_water_molecules(x, y, z)
    inp = tempfile.NamedTemporaryFile(mode="w+", suffix="inp")
    inp.write("\n".join([
        f"tolerance {tol}",
        "filetype pdb",
        "output solvated.pdb\n",
        f"structure {directory}/H2O.pdb",
        f"  number {num_waters}",
        f"  inside box 0. 0. 0. {x} {y} {z}\n"
        "end structure\n",
    ]))
    for ion, concentration in kwargs.items():
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
    return inp


if __name__ == "__main__":
    import subprocess
    import argparse

    ap = argparse.ArgumentParser("Solution", description="Disolve molecules in water.")
    size_ap = ap.add_mutually_exclusive_group(required=True)
    soluble_ap = ap.add_mutually_exclusive_group(required=True)

    ap.add_argument("-d", "--soluble-dir", dest="directory", default="./mixture/", help="Looking for solubles in this directory.")
    ap.add_argument("-t", "--tolerance", type=float, default=2.0, help="Distance tolerance.")
    ap.add_argument("-o", "--out", default="./solvated.pdb", help="Output pdb file. Default is ./solvated.pdb")
    size_ap.add_argument("-s", "--size", type=float, nargs=3)
    soluble_ap.add_argument("-a", "--add", dest="solubles", nargs=2, action="append", help="Specify type and concentration.")

    arguments = ap.parse_args()

    x, y, z = arguments.size
    solubles = {k: float(v) for k, v in arguments.solubles}
    print(solubles)

    packmol_input = solvate(
        x, y, z,
        arguments.tolerance,
        arguments.directory,
        **solubles
    )

    subprocess.call("packmol", stdin=packmol_input)
