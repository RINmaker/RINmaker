# RINmaker

## Clone this repo

```bash
git clone -b dev https://github.com/RINmaker/RINmaker.git
```

```bash
cd RINmaker
git submodule update --init --recursive
```

## Usage

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
