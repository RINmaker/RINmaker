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
    std::string x;
    std::string y;
    std::string z;
    std::string bfactor_ca;
    std::string secondary_structure;

    int degree = 0;
};

struct rin::edge::impl final
{
public:
    std::string source;
    std::string target;
    std::string source_atom;
    std::string target_atom;
    std::string distance;
    std::string energy;
    std::string angle;
    std::string interaction;
    std::string orientation;
    std::string donor;
    std::string cation;
    std::string positive;
};