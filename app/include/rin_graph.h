#pragma once

#pragma warning(push, 0)

#include <pugixml.hpp>

#pragma warning(pop)

#include <string>
#include <vector>
#include <memory>
#include <filesystem>
#include <unordered_map>

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
    std::string const& get_source_id() const;

    [[nodiscard]]
    std::string const& get_target_id() const;

    [[nodiscard]]
    std::string const& get_distance() const;

    [[nodiscard]]
    std::string const& get_energy() const;

    [[nodiscard]]
    std::string const& get_interaction() const;

    [[nodiscard]]
    std::string const& get_source_atom() const;

    [[nodiscard]]
    std::string const& get_target_atom() const;

    [[nodiscard]]
    std::string const& get_angle() const;

    [[nodiscard]]
    std::string const& get_donor() const;

    [[nodiscard]]
    std::string const& get_cation() const;

    [[nodiscard]]
    std::string const& get_positive() const;

    [[nodiscard]]
    std::string const& get_orientation() const;

    void append_to(pugi::xml_node& rin, bool with_metadata = false) const;
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

    node& operator++();

    [[nodiscard]]
    int get_degree() const;

    [[nodiscard]]
    std::string const& get_id() const;

    void append_to(pugi::xml_node& graphml, bool with_metadata = false) const;
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
        std::vector<chemical_entity::aminoacid> const& aminoacids,
        std::vector<std::shared_ptr<bond::base const>> const& bonds);

    graph(graph const& other);

    ~graph();

    void write_to_file(std::filesystem::path const& out_path) const;

    [[nodiscard]]
    std::string get_name() const;

    [[nodiscard]]
    std::vector<edge> const& get_edges() const;

    [[nodiscard]]
    std::unordered_map<std::string, node> const& get_nodes() const;
};
}