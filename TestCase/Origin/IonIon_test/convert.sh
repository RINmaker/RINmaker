#!/bin/bash
 
in_type="cif"
out_type="pdb"
read filename


obabel -i "${in_type}" "${filename}.cif" -o "${out_type}" -O "${filename}_pdb"
