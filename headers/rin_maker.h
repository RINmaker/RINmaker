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

#include "chemical_entities.h"
#include "bond_queries.h"
#include "bond_network.h"

#include "utils/spatial/kdtree.h"

namespace rin_maker {

class base {
protected:
    std::vector<entities::aminoacid*> _aminoacids;
    network _rin_network;
    rin::graph _rin_graph;

    explicit base(std::filesystem::path const& pdb_path);

public:
    virtual ~base();
};

class all_bonds final : public base {
private:
    kdtree<entities::atom, 3> _hdonor_tree, _vdw_tree;
    std::vector<entities::atom const*> _hacceptor_vector, _vdw_vector, _cation_vector;

    kdtree<entities::ring, 3> _ring_tree, _pication_ring_tree;
    std::vector<entities::ring const*> _ring_vector, _pication_ring_vector;

    kdtree<entities::ionic_group, 3> _positive_ion_tree;
    std::vector<entities::ionic_group const*> _negative_ion_vector;

public:
    explicit all_bonds(std::filesystem::path const& pdb_path);
};

class carbon : public base {
protected:
    kdtree<entities::atom, 3> _carbon_tree;
    std::vector<entities::atom const*> _carbon_vector;

    explicit carbon(std::filesystem::path const& pdb_path, std::function<entities::atom const*(entities::aminoacid const&)> const& getter);
};

class alpha_carbon final : public carbon {
public:
    explicit alpha_carbon(std::filesystem::path const& pdb_path)
            : carbon(
            pdb_path, [](entities::aminoacid const& res) -> entities::atom const* { return res.ca(); }) {
        log_manager::main()->info("alpha carbons: {}", _carbon_vector.size());
    }
};

class beta_carbon final : public carbon {
public:
    explicit beta_carbon(std::filesystem::path const& pdb_path)
            : carbon(pdb_path, [](entities::aminoacid const& res) -> entities::atom const* { return res.cb(); }) {
        log_manager::main()->info("beta carbons: {}", _carbon_vector.size());
    }
};
}
