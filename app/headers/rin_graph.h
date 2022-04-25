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

class generic_bond;
}

namespace rin
{

using std::string;

class edge
{
private:
    struct impl;
    impl* pimpl;

public:
    explicit edge(bond::ss const& bond);

    explicit edge(bond::vdw const& bond);

    explicit edge(bond::ionic const& bond);

    explicit edge(bond::hydrogen const& bond);

    explicit edge(bond::pication const& bond);

    explicit edge(bond::pipistack const& bond);

    explicit edge(bond::generic_bond const& bond);

    edge(edge const&);

    ~edge();

    [[nodiscard]]
    std::string const& source_id() const;

    [[nodiscard]]
    std::string const& target_id() const;

    [[nodiscard]]
    std::string const& distance() const;

    [[nodiscard]]
    std::string const& energy() const;

    [[nodiscard]]
    std::string const& interaction() const;

    [[nodiscard]]
    std::string const& source_atom() const;

    [[nodiscard]]
    std::string const& target_atom() const;

    [[nodiscard]]
    std::string const& angle() const;

    [[nodiscard]]
    std::string const& donor() const;

    [[nodiscard]]
    std::string const& cation() const;

    [[nodiscard]]
    std::string const& positive() const;

    [[nodiscard]]
    std::string const& orientation() const;

    void append_to(pugi::xml_node& rin, bool metadata);
};

class node
{
private:
    struct impl;
    impl* pimpl;

public:
    explicit node(chemical_entity::aminoacid const& res);

    node(node const& other);

    ~node();

    void inc_degree();

    [[nodiscard]]
    int degree() const;

    [[nodiscard]]
    string const& get_id() const;

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
    graph(string name, parameters const& params, vector<aminoacid const*> const& aminoacids,
          vector<std::shared_ptr<bond::base const>> const& bonds);

    void write_to_file(fs::path const& out_path) const;

    string name() const
    { return _name; }

    std::vector<edge> get_edges() const
    { return _edges; }

    std::unordered_map<std::string, node> get_nodes() const
    {
        std::unordered_map<std::string, node> out;
        for (const auto& i: _nodes)
            out.insert_or_assign(i.first, i.second);
        return out;
    }
};
}