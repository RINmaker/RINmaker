#pragma once

#include <string>
#include <vector>
#include <list>

#include <map>
#include <unordered_map>

#include <functional>
#include <filesystem>
#include <fstream>
#include <exception>
#include <memory>

#include "chemical_entity.h"
#include "bond_queries.h"
#include "rin_network.h"

#include "spatial/kdtree.h"

namespace rin_maker
{

namespace fs = std::filesystem;

using lm = log_manager;

using std::vector;
using std::list;
using std::string;
using std::unordered_map;
using std::map;
using std::function;

using chemical_entity::aminoacid;
using chemical_entity::atom;
using chemical_entity::ring;
using chemical_entity::ionic_group;

using std::unique_ptr;

struct rin_maker final
{
private:
    vector<aminoacid*> _aminoacids;

    // stuff for non-convalent bonds
    kdtree<atom, 3> _hdonor_tree, _vdw_tree;
    vector<atom const*> _hacceptor_vector, _vdw_vector, _cation_vector;

    kdtree<ring, 3> _ring_tree, _pication_ring_tree;
    vector<ring const*> _ring_vector, _pication_ring_vector;

    kdtree<ionic_group, 3> _positive_ion_tree;
    vector<ionic_group const*> _negative_ion_vector;

    // stuff for backbone
    kdtree<atom, 3> _alpha_carbon_tree, _beta_carbon_tree;
    vector<atom const*> _alpha_carbon_vector, _beta_carbon_vector;

public:
    explicit rin_maker(fs::path const& pdb_path);

    ~rin_maker()
    { for (auto* res: _aminoacids) delete res; }

    rin::graph operator()(parameters::interaction_type interaction_type, parameters::policy network_policy) const;
};
}