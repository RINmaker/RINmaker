#pragma once

#include "rin_maker.h"

#include <memory>
#include <vector>
#include <string>

#include "chemical_entity.h"

struct rin::maker::impl
{
public:
    std::vector<chemical_entity::aminoacid> aminoacids;

    kdtree<chemical_entity::atom, 3> hdonor_tree, vdw_tree;
    std::vector<chemical_entity::atom> hacceptor_vector, vdw_vector, cation_vector;

    kdtree<chemical_entity::ring, 3> ring_tree, pication_ring_tree;
    std::vector<chemical_entity::ring> ring_vector, pication_ring_vector;

    kdtree<chemical_entity::ionic_group, 3> positive_ion_tree;
    std::vector<chemical_entity::ionic_group> negative_ion_vector;

    kdtree<chemical_entity::atom, 3> alpha_carbon_tree, beta_carbon_tree;
    std::vector<chemical_entity::atom> alpha_carbon_vector, beta_carbon_vector;

    // ss bonds are directly parsed, not computed by us
    std::vector<std::shared_ptr<bond::ss const>> ss_bonds;

    std::string pdb_name;
};
