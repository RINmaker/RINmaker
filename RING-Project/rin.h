#pragma once

#include <string>
#include <array>
#include <queue>
#include <unordered_map>
#include <pugixml/pugixml.hpp>

// entity has to #include "rin"
// #include "entity.h"
namespace entities { class aminoacid; }

namespace bonds
{
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
class edge
{
private:
    std::string _source, _target;
    std::string _distance, _energy;

    std::string _interaction;
    std::string _source_atom, _target_atom;

    std::string _angle;
    std::string _donor;
    std::string _cation;
    std::string _positive;
    std::string _orientation;

public:
    edge(bonds::ss const& bond);
    edge(bonds::vdw const& bond);
    edge(bonds::ionic const& bond);
    edge(bonds::hydrogen const& bond);
    edge(bonds::pication const& bond);
    edge(bonds::pipistack const& bond);
    edge(bonds::generico const& bond);

public:
    std::string const& source_id() const { return _source; }
    std::string const& target_id() const { return _target; }

    std::string const& distance() const { return _distance; }
    std::string const& energy() const { return _energy; }

    std::string const& interaction() const { return _interaction; }
    std::string const& source_atom() const { return _source_atom; }
    std::string const& target_atom() const { return _target_atom; }

    std::string const& angle() const { return _angle; }
    std::string const& donor() const { return _donor; }
    std::string const& cation() const { return _cation; }
    std::string const& positive() const { return _positive; }
    std::string const& orientation() const { return _orientation; }

    void append_to(pugi::xml_node& rin, bool metadata);
};

class node
{
private:
    std::string _id;
    std::string _pdb_name;
    std::string _chain;
    std::string _seq;
    std::string _name;
    std::string _x, _y, _z;
    std::string _bfactor;
    std::string _secondary;

private:
    int _degree = 0;

public:
    explicit node(entities::aminoacid const& res);

    std::string const& get_id() const { return _id; }
    int& degree() { return _degree; }

    void append_to(pugi::xml_node& graphml, bool metadata) const;
};

class graph
{
private:
    std::unordered_map<std::string, node> nodes;
    std::queue<edge> edges;

    edge pop_edge();

public:
    void push(edge const& e) { edges.push(e); }
    void insert(node const& n) { nodes.insert({ n.get_id(), n }); }

public:
    void consume_to_xml();
    std::vector<edge> get_edges();
    std::unordered_map<std::string, node> get_nodes();
};
}