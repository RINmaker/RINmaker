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

namespace rin_maker {

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

////////////////////////////////////////////////////////////////////////

class base {
protected:
    vector<aminoacid*> _aminoacids;
    network _rin_network;

    explicit base(fs::path const& pdb_path);

public:
    virtual ~base() {
        for (auto* res: _aminoacids)
            delete res;
    }

    rin::graph get_graph(parameters::interaction_type inter) const;
};

class all_bonds final : public base {
private:
    double _query_dist_hbond, _query_dist_vdw, _query_dist_ionic, _query_dist_pipi, _query_dist_pica;

    kdtree<atom, 3> _hdonor_tree, _vdw_tree;
    vector<atom const*> _hacceptor_vector, _vdw_vector, _cation_vector;

    kdtree<ring, 3> _ring_tree, _pication_ring_tree;
    vector<ring const*> _ring_vector, _pication_ring_vector;

    kdtree<ionic_group, 3> _positive_ion_tree;
    vector<ionic_group const*> _negative_ion_vector;

public:
    explicit all_bonds(fs::path const& pdb_path);
};

typedef function<atom const*(aminoacid const&)> const& carbon_getter;

class backbone : public base {
private:
    double _query_dist_carbon;

protected:
    kdtree<atom, 3> _carbon_tree;
    vector<chemical_entity::atom const*> _carbon_vector;

    explicit backbone(fs::path const& pdb_path, carbon_getter getter);
};

class alpha_carbon final : public backbone {
public:
    explicit alpha_carbon(fs::path const& pdb_path)
            : backbone(
            pdb_path, [](aminoacid const& res) -> atom const* {
                return res.ca();
            }) { lm::main()->info("alpha carbons: {}", _carbon_vector.size()); }
};

class beta_carbon final : public backbone {
public:
    explicit beta_carbon(fs::path const& pdb_path)
            : backbone(
            pdb_path, [](aminoacid const& res) -> atom const* {
                return res.cb();
            }) { lm::main()->info("beta carbons: {}", _carbon_vector.size()); }
};

unique_ptr<base const> make_instance(parameters::policy policy, fs::path const& path);
}
