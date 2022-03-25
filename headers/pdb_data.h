#pragma once

#include <vector>
#include <memory>

#include "runtime_params.h"
#include "bond_network.h"

#include "spatial/kdtree.h"
#include "graphml_output.h"

#include "chemical_entity.h"
#include "log_manager.h"

class pdb_data {
private:
    std::vector<chemical_entity::aminoacid *> _aminoacids;

    std::vector<chemical_entity::atom const *> _hacceptors, _vdws, _cas, _cbs, _cations;
    std::vector<chemical_entity::ring const *> _rings;
    std::vector<chemical_entity::ring const *> _pication_rings;
    std::vector<chemical_entity::ionic_group const *> _negatives;

private:
    kdtree<chemical_entity::atom, 3> _hdonors_tree, _vdws_tree, _cas_tree, _cbs_tree;
    kdtree<chemical_entity::ring, 3> _rings_tree;
    kdtree<chemical_entity::ring, 3> _pication_rings_tree;
    kdtree<chemical_entity::ionic_group, 3> _positives_tree;

private:
    network net;
    rin::graph g;

private:
    // parse atoms and other chemical stuff from pdb file.
    void parse_file();

public:
    rin::graph &get_graph() { return g; }

public:
    explicit pdb_data(bool exportOutput = true);

    ~pdb_data() {
        for (auto *res: _aminoacids)
            delete res;
    }
};
