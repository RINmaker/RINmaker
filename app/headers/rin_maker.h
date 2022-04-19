#pragma once

#include <vector>
#include <utility>
#include <filesystem>

#include "rin_graph.h"
#include "chemical_entity.h"

#include "log_manager.h"
#include "spatial/kdtree.h"

namespace rin
{
namespace fs = std::filesystem;

using lm = log_manager;

using std::vector;

using chemical_entity::aminoacid;
using chemical_entity::atom;
using chemical_entity::ring;
using chemical_entity::ionic_group;

struct parameters;

typedef std::pair<uint32_t, std::string> numbered_line_t;
struct maker final
{
private:
    vector<aminoacid const*> _aminoacids;

    kdtree<atom, 3> _hdonor_tree, _vdw_tree;
    vector<atom const*> _hacceptor_vector, _vdw_vector, _cation_vector;

    kdtree<ring, 3> _ring_tree, _pication_ring_tree;
    vector<ring const*> _ring_vector, _pication_ring_vector;

    kdtree<ionic_group, 3> _positive_ion_tree;
    vector<ionic_group const*> _negative_ion_vector;

    kdtree<atom, 3> _alpha_carbon_tree, _beta_carbon_tree;
    vector<atom const*> _alpha_carbon_vector, _beta_carbon_vector;

    // ss bonds are directly parsed, not computed by us
    vector<std::shared_ptr<bond::ss const>> _ss_bonds;

    std::string _pdb_name;

public:
    explicit maker(std::string const& pdb_name, std::vector<numbered_line_t>::iterator begin, std::vector<numbered_line_t>::iterator end);
    ~maker();

    rin::graph operator()(parameters const& params) const;
};
}