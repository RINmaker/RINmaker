# RINmaker

## CLI Usage

Assuming `RINmaker` is on `$PATH`:

```bash
RINmaker [OPTIONS] SUBCOMMAND
```

A brief help is available with `-h` or `--help`.

### Options

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

### Subcommands

#### `rin` options

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

#### `cmap` options

|    Param.    | Default | Meaning                                                                                                                              |
|:------------:|:-------:|--------------------------------------------------------------------------------------------------------------------------------------|
|   `--type`   |  `ca`   | <ul><li>`ca`: use alpha carbons.</li><li>`cb`: use beta carbons.</li></ul>                                                           |       
| `--distance` |   3.5   | Query distance between alpha/beta carbons.                                                                                           |


### Example

```bash
RINmaker -i 6j8j.pdb -o testrun.graphml rin
```

Will parse the first model in the *6j8j* pdb. 

## Clone this repo

```bash
git clone -b dev https://github.com/RINmaker/RINmaker.git
```

```bash
cd RINmaker
git submodule update --init --recursive
```

```bash
mkdir ~/.RINmaker
tar -xf monomers.tar.gz -C ~/.RINmaker/monomers
```


## Configure project:

```bash
cmake -S . -B build
```

## App

#### Build

```bash
cmake --build build --target RINmaker
```

#### Run

```bash
./build/app/RINmaker -h
```

## Test

#### Build

```bash
cmake --build build --target RINmaker_test
```

#### Run

```bash
./build/test/RINmaker_test
```
