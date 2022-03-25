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

namespace che = chemical_entity;
namespace fs = std::filesystem;

using lm = log_manager;

using std::vector;
using std::function;

class base {
protected:
    vector<che::aminoacid*> _aminoacids;
    network _rin_network;

    explicit base(fs::path const& pdb_path);

public:
    virtual ~base() {
        for (auto* res: _aminoacids)
            delete res;
    }

    rin::graph get_graph() const;
};

class all_bonds final : public base {
private:
    double _query_dist_hbond, _query_dist_vdw, _query_dist_ionic, _query_dist_pipi, _query_dist_pica;

    kdtree<che::atom, 3> _hdonor_tree, _vdw_tree;
    vector<che::atom const*> _hacceptor_vector, _vdw_vector, _cation_vector;

    kdtree<che::ring, 3> _ring_tree, _pication_ring_tree;
    vector<che::ring const*> _ring_vector, _pication_ring_vector;

    kdtree<che::ionic_group, 3> _positive_ion_tree;
    vector<che::ionic_group const*> _negative_ion_vector;

public:
    explicit all_bonds(fs::path const& pdb_path);
};

typedef function<che::atom const*(che::aminoacid const&)> const& carbon_getter;

class backbone : public base {
private:
    double _query_dist_carbon;

protected:
    kdtree<che::atom, 3> _carbon_tree;
    vector<chemical_entity::atom const*> _carbon_vector;

    explicit backbone(fs::path const& pdb_path, carbon_getter getter);
};

class alpha_carbon final : public backbone {
public:
    explicit alpha_carbon(fs::path const& pdb_path)
            : backbone(
            pdb_path, [](che::aminoacid const& res) -> che::atom const* {
                return res.ca();
            }) { lm::main()->info("alpha carbons: {}", _carbon_vector.size()); }
};

class beta_carbon final : public backbone {
public:
    explicit beta_carbon(fs::path const& pdb_path)
            : backbone(
            pdb_path, [](che::aminoacid const& res) -> che::atom const* {
                return res.cb();
            }) { lm::main()->info("beta carbons: {}", _carbon_vector.size()); }
};

base const* build(parameters::policy policy, fs::path const& path);
}
