#pragma once

#pragma warning(push, 0)

#include <gemmi/pdb.hpp>

#pragma warning(pop)

#include <memory>

#include "rin_graph.h"

namespace rin
{
struct parameters;

struct maker final
{
private:
    struct impl;
    std::shared_ptr<impl const> pimpl;

public:
    maker(gemmi::Model const& model, gemmi::Structure const& structure);

    ~maker();

    rin::graph operator()(parameters const& params) const;
};
}