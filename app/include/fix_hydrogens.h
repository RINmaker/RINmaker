#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <iostream>

#pragma warning(push, 0)

#include <gemmi/mmcif.hpp>     // for make_structure
#include <gemmi/polyheur.hpp>  // for setup_entities
#include <gemmi/modify.hpp>    // for remove_hydrogens
#include <gemmi/to_cif.hpp>    // for write_cif_to_file
#include <gemmi/to_mmcif.hpp>  // for make_mmcif_document
#include <gemmi/to_pdb.hpp>    // for write_pdb
#include <gemmi/monlib.hpp>    // for MonLib, read_monomer_lib
#include <gemmi/topo.hpp>      // for Topo
#include <gemmi/fstream.hpp>   // for Ofstream
#include <gemmi/placeh.hpp>    // for prepare_topology
#include <gemmi/read_coor.hpp> // for read_structure_gz
#include <gemmi/cif.hpp>

#pragma warning(pop)

void fix_hydrogens(gemmi::Structure& structure, gemmi::HydrogenChange what);