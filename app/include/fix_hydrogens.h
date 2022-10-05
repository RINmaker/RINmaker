#pragma warning(push, 0)

#include <gemmi/polyheur.hpp>   // for setup_entities
#include <gemmi/modify.hpp>     // for remove_hydrogens
#include <gemmi/monlib.hpp>     // for MonLib, read_monomer_lib
#include <gemmi/topo.hpp>       // for Topo
#include <gemmi/placeh.hpp>     // for prepare_topology
#include <gemmi/cif.hpp>        // for read_file

#pragma warning(pop)

void fix_hydrogens(gemmi::Structure& structure, gemmi::HydrogenChange what);