#pragma once

#include <string>
#include <vector>
#include <list>
#include <array>
#include <queue>
#include <unordered_map>
#include <filesystem>
#include <pugixml.hpp>

#include "rin_params.h"

namespace chemical_entity
{
class aminoacid;
}

namespace bond
{
class base;

class hydrogen;

class ss;

class vdw;

class pication;

class pipistack;

class ionic;

class generico;
}

namespace rin
{

using std::string;

class edge
{
private:
    string _source, _target, _source_atom, _target_atom;
    string _distance, _energy, _angle, _interaction, _orientation;
    string _donor, _cation, _positive;

public:
    explicit edge(bond::ss const& bond);

    explicit edge(bond::vdw const& bond);

    explicit edge(bond::ionic const& bond);

    explicit edge(bond::hydrogen const& bond);

    explicit edge(bond::pication const& bond);

    explicit edge(bond::pipistack const& bond);

    explicit edge(bond::generico const& bond);

public:
    [[nodiscard]]
    string const& source_id() const
    { return _source; }

    [[nodiscard]]
    string const& target_id() const
    { return _target; }

    [[nodiscard]]
    string const& distance() const
    { return _distance; }

    [[nodiscard]]
    string const& energy() const
    { return _energy; }

    [[nodiscard]]
    string const& interaction() const
    { return _interaction; }

    [[nodiscard]]
    string const& source_atom() const
    { return _source_atom; }

    [[nodiscard]]
    string const& target_atom() const
    { return _target_atom; }

    [[nodiscard]]
    string const& angle() const
    { return _angle; }

    [[nodiscard]]
    string const& donor() const
    { return _donor; }

    [[nodiscard]]
    string const& cation() const
    { return _cation; }

    [[nodiscard]]
    string const& positive() const
    { return _positive; }

    [[nodiscard]]
    string const& orientation() const
    { return _orientation; }

    void append_to(pugi::xml_node& rin, bool metadata);
};

class node
{
private:
    string _id;
    string _pdb_name;
    string _chain;
    string _seq;
    string _name;
    string _x, _y, _z;
    string _bfactor;
    string _secondary;

    int _degree = 0;

public:
    explicit node(chemical_entity::aminoacid const& res);

    void inc_degree()
    { ++_degree; }

    [[nodiscard]]
    int degree() const
    { return _degree; }

    [[nodiscard]]
    string const& get_id() const
    { return _id; }

    void append_to(pugi::xml_node& graphml, bool with_metadata) const;
};

using std::unordered_map;
using std::queue;
using std::vector;
using std::list;
using chemical_entity::aminoacid;
namespace fs = std::filesystem;

class graph
{
private:
    string _name;
    rin::parameters _params;
    unordered_map<string, node> _nodes;
    vector<edge> _edges;

public:
    graph(string name, parameters const& params, vector<aminoacid const*> const& aminoacids, vector<std::shared_ptr<bond::base const>> const& bonds);

    void write_to_file(fs::path const& out_path) const;

    string name() const { return _name; }

    std::vector<edge> get_edges() const { return _edges; }

    std::unordered_map<std::string, node> get_nodes() const
    {
        std::unordered_map<std::string, node> out;
        for (const auto& i: _nodes)
            out.insert_or_assign(i.first, i.second);
        return out;
    }
};
}