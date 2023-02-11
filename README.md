# RINmaker

## Table of Contents

* [CLI Usage](#usage)
    * [Example](#example)
    * [Help](#help)
    * [Options](#options)
    * [Subcommands](#subcommands)
        * [rin](#rin)
        * [cmap](#cmap)
* [Build instructions](#build)
    * [App](#app)
    * [Tests](#tests)

## CLI Usage <a name="usage"></a>

Assuming `RINmaker` is on `$PATH`:

```bash
RINmaker [OPTIONS] SUBCOMMAND
```

A brief help is available with `-h` or `--help`.

### Example <a name="example"></a>

```bash
RINmaker -i 6j8j.pdb -o testrun.graphml rin
```

Will parse the first model in the *6j8j* pdb.

### Help <a name="help"></a>

```
RINmaker v0.1.3 build Jan 10 2023 17:29:50 (Linux) 
(C) 2020-23 Ca' Foscari University of Venice

Usage: ./RINmaker [OPTIONS] SUBCOMMAND

Options:
  -h,--help                                                     Print this help message and exit
  -H,--help-expanded                                            Print this help message (expanded) and exit
  -i,--input TEXT:FILE REQUIRED                                 Path to .pdb or .cif file
  -o,--output TEXT REQUIRED                                     Output file (or directory if -d flag is specified)
  -d                                                            Use -o argument as a directory
  -l,--log TEXT=./main.txt                                      Log file
  -v,--verbose                                                  Log also to stdout
  -n,--no-hydrogen                                              Skip hydrogen fixing
  -w,--keep-water                                               Keep water residues
  -s,--sequence-separation INT:POSITIVE=3                       Minimum sequence separation
  --illformed ENUM:{fail,kall,kres,sres}=sres                   Behaviour in case of malformed ring or ionic group

Subcommands:
rin
  Compute the residue interaction network
  Options:
    --policy ENUM:{all,multiple,one}=all                          Affects which edges are kept per pair of aminoacids
    --hydrogen-bond FLOAT:POSITIVE=3.5                            Query distance for hydrogen bonds
    --vdw-bond FLOAT=0.5                                          Surface distance for vdw bonds
    --ionic-bond FLOAT:POSITIVE=4                                 Query distance for ionic bonds
    --pication-bond FLOAT:POSITIVE=5                              Query distance for cation-pi bonds
    --pipistack-bond FLOAT:POSITIVE=6.5                           Query distance for pi-pi stackings
    --h-bond-realistic                                            Keep only MC-MC hydrogen bonds with minimum energy
    --h-bond-angle FLOAT:POSITIVE=63                              Angle for hydrogen bonds
    --pication-angle FLOAT:POSITIVE=45                            Angle for cation-pi bonds
    --pipistack-normal-normal FLOAT:POSITIVE=30                   Angle range from normal to normal for pi-pi stackings
    --pipistack-normal-centre FLOAT:POSITIVE=60                   Angle range from normal to centre for pi-pi stackings

cmap
  Compute the contact map of the protein
  Options:
    --type ENUM:{ca,cb}=ca                                        Type of contact map (alpha/beta carbon)
    --distance FLOAT:POSITIVE=6                                   Query distance between alpha/beta carbons
```

### Options <a name="options"></a>

|         Param.          | Abbr. |     Default     | Meaning                                                                                                                                                                                                                    |
|:-----------------------:|:-----:|:---------------:|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|        `--help`         | `-h`  |      `all`      | <ul><li>`all`: keep all bonds.</li><li>`one`: for each pair of residues, output the least energetic bond.</li><li>`multiple`: for each pair of residues, for each type of bond, output the least energetic bond.</li></ul> |       
|    `--help-expanded`    | `-H`  |       3.5       | Maximum distance between atoms (in ångström) for a pair donor-acceptor to be tested for hydrogen bond.                                                                                                                     |       
|        `--input`        | `-i`  | <u>REQUIRED</u> | Maximum distance (in ångström) between the _surfaces_ of the spheres marked by two atoms' Van der Waals radii.                                                                                                             |
|       `--output`        | `-o`  | <u>REQUIRED</u> | If `-d` is not used, then it represents the path to the output file. If `-d` is used, it will be used as a directory for the resulting graphs.                                                                             |
|                         | `-d`  |     not set     | It's a flag. By default, it is not set and the program will analyse only the first model of the protein. If used, then all the models in the PDB/mmCIF will be analysed instead and multiple graphs will be produced.      |
|         `--log`         | `-l`  |  "./main.txt"   | Maximum distance between two aromatic rings (in ångström; centre of mass is considered) to be tested for bonding conditions.                                                                                               |
|       `--verbose`       | `-v`  |     not set     | It's a flag. If used, logs will be mirrored to standard output.                                                                                                                                                            |
|     `--no-hydrogen`     | `-n`  |     not set     | Skip hydrogen fixing.                                                                                                                                                                                                      |
|     `--keep-water`      | `-w`  |     not set     | Keep water residues                                                                                                                                                                                                        |
| `--sequence-separation` | `-s`  |        3        | Minimum sequence separation                                                                                                                                                                                                |
|      `--illformed`      | `-f`  |     `sres`      | <ul><li>`kall`: keep everything.</li><li>`kres`: keep the residue _without_ considering the malformed part.</li><li>`sres`: skip the residue altogether.</li><li>`fail`: halt with error.</li></ul>                        |

### Subcommands <a name="subcommands"></a>

#### `rin` options <a name="rin"></a>

|           Param.            | Default | Meaning                                                                                                                                                                                                                    |
|:---------------------------:|:-------:|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|         `--policy`          |  `all`  | <ul><li>`all`: keep all bonds.</li><li>`one`: for each pair of residues, output the least energetic bond.</li><li>`multiple`: for each pair of residues, for each type of bond, output the least energetic bond.</li></ul> |       
|      `--hydrogen-bond`      |   3.5   | Maximum distance between atoms (in ångström) for a pair donor-acceptor to be tested for hydrogen bond.                                                                                                                     |       
|        `--vdw-bond`         |   0.5   | Maximum distance (in ångström) between the _surfaces_ of the spheres marked by two atoms' Van der Waals radii.                                                                                                             |
|       `--ionic-bond`        |   4.0   | Maximum distance between two ionic groups (in ångström; centre of mass is considered) to be tested for bonding conditions.                                                                                                 |
|      `--pication-bond`      |   5.0   | Maximum distance between a cation and an aromatic ring (in ångström; centre of mass is considered) to be considered bond candidates.                                                                                       |
|     `--pipistack-bond`      |   6.5   | Maximum distance between two aromatic rings (in ångström; centre of mass is considered) to be tested for bonding conditions.                                                                                               |
|    `--h-bond-realistic`     | not set | It's a flag; if used, it keeps only MC-MC hydrogen bonds with minimum energy _before_ filtering by policy.                                                                                                                 |
|      `--h-bond-angle`       |   63    | Angle (degrees) for hydrogen bonds.                                                                                                                                                                                        |
|     `--pication-angle`      |   45    | Angle for cation-pi bonds.                                                                                                                                                                                                 |
| `--pipistack-normal-normal` |   30    | Angle range from normal to normal for pi-pi stackings.                                                                                                                                                                     |
| `--pipistack-normal-centre` |   60    | Angle range from normal to centre for pi-pi stackings.                                                                                                                                                                     |

#### `cmap` options <a name="cmap"></a>

|    Param.    | Default | Meaning                                                                                                                              |
|:------------:|:-------:|--------------------------------------------------------------------------------------------------------------------------------------|
|   `--type`   |  `ca`   | <ul><li>`ca`: use alpha carbons.</li><li>`cb`: use beta carbons.</li></ul>                                                           |       
| `--distance` |   3.5   | Query distance between alpha/beta carbons.                                                                                           |


## Build instructions <a name="build"></a>

You will need: 
* `cmake` 
* A suitable C++ compiler. For *nix systems `g++`; for Windows, you can use `MSVC`.

Clone and initialize the project:

```bash
git clone -b dev https://github.com/RINmaker/RINmaker.git
cd RINmaker
git submodule update --init --recursive
cmake -S . -B build
```
### App <a name="app"></a>

```bash
cmake --build build --target RINmaker
```
The application's executable will be located at: `./build/app/RINmaker`.

### Tests <a name="tests"></a>

```bash
cmake --build build --target RINmaker_test
```

To run the _full_ test suite (assuming you are in the project's root):

```bash
./build/test/RINmaker_test
```
