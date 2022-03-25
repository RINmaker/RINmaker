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

#include "chemical_entity.h"
#include "bond_queries.h"
#include "bond_network.h"

#include "spatial/kdtree.h"

namespace rin_maker {

class base {
protected:
    std::vector<chemical_entity::aminoacid*> _aminoacids;
    network _rin_network;
    rin::graph _rin_graph;

    explicit base(std::filesystem::path const& pdb_path);

public:
    virtual ~base();
};

class all_bonds final : public base {
private:
    kdtree<chemical_entity::atom, 3> _hdonor_tree, _vdw_tree;
    std::vector<chemical_entity::atom const*> _hacceptor_vector, _vdw_vector, _cation_vector;

    kdtree<chemical_entity::ring, 3> _ring_tree, _pication_ring_tree;
    std::vector<chemical_entity::ring const*> _ring_vector, _pication_ring_vector;

    kdtree<chemical_entity::ionic_group, 3> _positive_ion_tree;
    std::vector<chemical_entity::ionic_group const*> _negative_ion_vector;

public:
    explicit all_bonds(std::filesystem::path const& pdb_path);
};

typedef std::function<chemical_entity::atom const*(chemical_entity::aminoacid const&)> const& carbon_getter;
class backbone : public base {
protected:
    kdtree<chemical_entity::atom, 3> _carbon_tree;
    std::vector<chemical_entity::atom const*> _carbon_vector;

    explicit backbone(std::filesystem::path const& pdb_path, carbon_getter getter);
};

class alpha_carbon final : public backbone {
public:
    explicit alpha_carbon(std::filesystem::path const& pdb_path)
            : backbone(
            pdb_path, [](chemical_entity::aminoacid const& res) -> chemical_entity::atom const* {
                return res.ca();
            }) { log_manager::main()->info("alpha carbons: {}", _carbon_vector.size()); }
};

class beta_carbon final : public backbone {
public:
    explicit beta_carbon(std::filesystem::path const& pdb_path)
            : backbone(pdb_path, [](chemical_entity::aminoacid const& res) -> chemical_entity::atom const* { return res.cb(); }) {
        log_manager::main()->info("beta carbons: {}", _carbon_vector.size());
    }
};
}
