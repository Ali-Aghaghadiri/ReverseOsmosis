**KEEP IN MIND THAT THIS REPOSITORY IS UNDER DEVELOPMENT**

# Dependency

Install dependencies using `pip`:
```shell
pip install -r requirements.txt
```
Also `packmol` and `VMD` are necessary.

> This project is powered by *ASE*, *VMD* and *packmol*.


# Crystal Structure

To generate the crystal of interest (in our case $Ti_3AlCN$) use `MyCrystal.py` script.
```shell
python3 MyCrystal.py --super-cell 16 16 1 --mx --mono --show cube
```
This command line creates a super cell of $Ti_3CN$ MXene phase by 16 units in $x$ and $y$ directions and 1 unit in $z$ direction. Then represents the crystal in cubic form.

We use this script to design *piston* and *membrane* parts of the system.
For creating the *membrane* we use `--drill` option:

```shell
python3 MyCrystal.py --super-cell 16 16 1 --mx --mono --show cube --drill 5.5
```
This will *drill* a hole in center of the crystal with radius of `5.5` $\mathring {\mathrm A}$. Or with the diameter of ~ 1 $nm$.

More details is accessible via `-h` or `--help`:

```
Titanium Aluminum Cyanide MAX and MXene crystal maker.

options:
  -h, --help            show this help message and exit

Lattice:
  Manipulate lattice parameters

  -a A                  Lattice parameter `a`
  -c C                  Lattice parameter `c`
  --super-cell SIZE SIZE SIZE
                        Repetition of unitcell

Phase:
  determine crystal phase.

  --max                 MAX phase
  --mx                  MX phase
  --mono-layer          Monolayer Structure

Output:
  Determine the output format.

  --format {xyz,pdb}    Output format
  -o DIRECTORY, --out DIRECTORY
                        Output directory
  --show {hex,cube}     View results

Manipulation:
  --drill DRILL         Make vacancy.
```

# Reverse Osmosis System

To create a RO system three main parts can be named:
- Membrane
- Piston
- Mixture

**Membrane** and **Piston** are created using `MyCrystal.py` and membrane has a hole in center.
For filling the chamber with desired solution we need to define molecules like $H_2O$ separately.
Our molecules are in `mixture` directory. Also the piston and membranes are in `piston` and `membrane` directories respectively.

## ROSystem.py

To create such system we can use `ROSystem.py` as followed:
```
python3 ROSystem.py -m membrane/Ti3CN.pdb -p piston/Ti3CN.pdb -s 130 -M Ca 3 -M Cl 3 -M Na 3 -o system/system.pdb --show
```

This command specifies membrane, piston and mixture compounds and relevant count for each type. And finally saves the result in a `pdb` file.
Then the script tries to make a **LAMMPS** `data` file from `pdb` file.
So there are two files generated with this command:
 - `system.pdb` *using packmol*
 - `system.data` *using VMD*

Alternatively we can create a RO system by just having a unit cell and desired size of the chamber and molecule concentrations:

```shell
 py ROSystem.py -c models/Ti3CN_mono.cif --size 50 50 100 -d 11 -M Ca 3 -M Cl 3 -M Na 3
```
This command will use the unit cell defined in a `cif` file and will create a $5\times5 nm$ piston and membrane with a hole of ~ $1 nm$. the chamber will be $10 nm$ long and will contain a solution of calcium sodium chlorid 3 Molar.

More details is accessible via `-h` or `--help`:
```
Reverse Osmosis System builder.

options:
  -h, --help            show this help message and exit
  -x EXECUTABLE, --executable EXECUTABLE
                        Packmol executable path.
  -o OUT, --out OUT     Output file.
  -t TOLERANCE, --tolerance TOLERANCE
                        Set tolerance.
  --show                View the packing result.
  -T TEMPERATURE, --temperature TEMPERATURE
                        Solution temperature.

Membrane:
  Specify characteristics of membrane.

  -m MEMBRANE, --membrane MEMBRANE
                        Membrane file.
  -d DRILL, --drill DRILL
                        Diameter of the membrane hole.

Piston:
  Specify characteristics of piston.

  -p PISTON, --piston PISTON
                        Piston file.

Mixture:
  Specify characteristics of mixture.

  -M MOLECULES MOLECULES, --molecule MOLECULES MOLECULES
                        Mixture compounds and concentration.

System:
  Specify characteristics of system.

  -s SEPARATION, --separation SEPARATION
                        Separation of the piston and membrane.
  -c CRYSTAL, --crystal CRYSTAL
                        Unit cell crystal file.
  --size SIZE SIZE SIZE
                        Size of the chamber.
```
