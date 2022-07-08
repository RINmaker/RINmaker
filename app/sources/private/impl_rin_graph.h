#pragma once

#include "rin_graph.h"

#include <string>
#include <vector>
#include <unordered_map>

#include "rin_params.h"

struct rin::graph::impl final
{
public:
    std::string name;
    rin::parameters params;
    std::unordered_map<std::string, rin::node> nodes;
    std::vector<rin::edge> edges;

    impl(std::string nm, rin::parameters const& pr) : name{std::move(nm)}, params{pr}
    {}
};

struct rin::node::impl final
{
public:
    std::string id;
    std::string pdb_name;
    std::string chain;
    std::string sequence_number;
    std::string name;
    std::string x, y, z;
    std::string bfactor_ca;
    std::string secondary_structure;

    int degree = 0;
};

struct rin::edge::impl final
{
public:
    std::string source, target, source_atom, target_atom;
    std::string distance, energy, angle, interaction, orientation;
    std::string donor, cation, positive;
};