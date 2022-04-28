#pragma once

#include <string>
#include <vector>
#include <memory>
#include <filesystem>
#include <unordered_map>
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
class edge
{
private:
    struct impl;
    std::shared_ptr<impl const> pimpl;

public:
    explicit edge(bond::ss const& bond);

    explicit edge(bond::vdw const& bond);

    explicit edge(bond::ionic const& bond);

    explicit edge(bond::hydrogen const& bond);

    explicit edge(bond::pication const& bond);

    explicit edge(bond::pipistack const& bond);

    explicit edge(bond::generic_bond const& bond);

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

    void append_to(pugi::xml_node& rin, bool with_metadata);
};

class node
{
private:
    struct impl;
    std::unique_ptr<impl> pimpl;

public:
    explicit node(chemical_entity::aminoacid const& res);

    node(node const& other);

    node& operator=(node const& rhs);

    ~node();

    void inc_degree();

    [[nodiscard]]
    int degree() const;

    [[nodiscard]]
    std::string const& get_id() const;

    void append_to(pugi::xml_node& graphml, bool with_metadata) const;
};

class graph
{
private:
    struct impl;
    std::shared_ptr<impl const> pimpl;

public:
    graph(
            std::string const& name,
            parameters const& params,
            std::vector<chemical_entity::aminoacid const*> const& aminoacids,
            std::vector<std::shared_ptr<bond::base const>> const& bonds);

    graph(graph const& other);

    ~graph();

    void write_to_file(std::filesystem::path const& out_path) const;

    std::string name() const;

    std::vector<edge> get_edges() const;

    std::unordered_map<std::string, node> get_nodes() const;
};
}