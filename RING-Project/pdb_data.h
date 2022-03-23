#pragma once

#include <vector>
#include <memory>

#include "parameters.h"
#include "network.h"

#include "kdtree.h"
#include "rin.h"

#include "entity.h"
#include "log_manager.h"

class pdb_data
{
private:
    std::shared_ptr<spdlog::logger> logger;

private:
    // amminoacidi della proteina
    std::vector<entities::aminoacid*> _aminoacids;

private:
    std::vector<entities::atom const*> _hacceptors, _vdws, _cas, _cbs, _cations;
    std::vector<entities::ring const*> _rings;
    std::vector<entities::ring const*> _pication_rings;
    std::vector<entities::ionic_group const*> _negatives;

private:
    kdtree<entities::atom, 3> _hdonors_tree, _vdws_tree, _cas_tree, _cbs_tree;
    kdtree<entities::ring, 3> _rings_tree;
    kdtree<entities::ring, 3> _pication_rings_tree;
    kdtree<entities::ionic_group, 3> _positives_tree;

private:
    network net;
    rin::graph g;

private:
    // parse atoms and other chemical stuff from pdb file.
    void parse_file();

public:
    rin::graph& get_graph()
    { return g; }

public:
    explicit pdb_data(bool exportOutput = true);

    ~pdb_data()
    {
        for (auto* res: _aminoacids)
            delete res;
    }
};
