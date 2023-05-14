#pragma once

#pragma warning(push, 0)

#include <cstdint> // apparently, gemmi does not include it properly (but it needs it).
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
    maker(gemmi::Model const& model, gemmi::Structure const& protein, rin::parameters const& params);

    ~maker();

    rin::graph operator()(parameters const& params) const;
};
}