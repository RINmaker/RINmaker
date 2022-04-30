#pragma once

#include <string>
#include <vector>
#include <array>
#include <memory>
#include <optional>

#include "config.h"

#include "rin_graph.h"
#include "ns_record.h"
#include "ns_secondary_structure.h"

#include "spatial/kdpoint.h"

namespace chemical_entity
{
class atom;

class ring;

class ionic_group;

class aminoacid : public kdpoint<3>
{
private:
    struct impl;
    std::shared_ptr<impl> pimpl;

public:
    class component
    {
    private:
        struct impl;
        std::shared_ptr<component::impl> pimpl;

    protected:
        explicit component(aminoacid const& res);

    public:
        [[nodiscard]]
        aminoacid res() const;
    };

    friend class component;

    aminoacid();

public:
    aminoacid(std::vector<record::atom> const& records, std::string const& pdb_name);

    ~aminoacid();

    [[nodiscard]]
    std::vector<atom> const& atoms() const;

    [[nodiscard]]
    std::string const& pdb_name() const;

    [[nodiscard]]
    std::optional<atom> const& ca() const;

    [[nodiscard]]
    std::optional<atom> const& cb() const;

    [[nodiscard]]
    std::optional<ring> const& primary_ring() const;

    [[nodiscard]]
    std::optional<ring> const& secondary_ring() const;

    [[nodiscard]]
    std::optional<ionic_group> const& positive_ionic_group() const;

    [[nodiscard]]
    std::optional<ionic_group> const& negative_ionic_group() const;

    [[nodiscard]]
    std::string const& name() const;

    [[nodiscard]]
    std::string const& chain_id() const;

    [[nodiscard]]
    std::string const& id() const;

    [[nodiscard]]
    int sequence_number() const;

    bool operator==(aminoacid const& rhs) const;

    bool operator!=(aminoacid const& rhs) const;

    [[nodiscard]]
    bool satisfies_minimum_separation(aminoacid const& aa, int minimum_separation = cfg::params::seq_sep) const;

    [[nodiscard]]
    explicit operator rin::node() const;

    void set_loop();

    void set_helix(record::helix const& record);

    void set_sheet(record::sheet_piece const& record);

    [[nodiscard]]
    std::string secondary_structure_id() const;
};

class atom final : public kdpoint<3>, public aminoacid::component
{
private:
    struct impl;
    std::shared_ptr<impl const> pimpl;

public:
    atom(record::atom const& record, aminoacid const& res);

    ~atom();

    [[nodiscard]]
    std::string const& name() const;

    [[nodiscard]]
    std::string const& symbol() const;

    [[nodiscard]]
    double temp_factor() const;

    [[nodiscard]]
    int charge() const;

    [[nodiscard]]
    bool is_in_a_positive_ionic_group() const;

    [[nodiscard]]
    bool is_in_a_negative_ionic_group() const;

    [[nodiscard]]
    bool is_a_hydrogen_donor() const;

    [[nodiscard]]
    int how_many_hydrogen_can_donate() const;

    [[nodiscard]]
    bool is_a_hydrogen_acceptor() const;

    [[nodiscard]]
    int how_many_hydrogen_can_accept() const;

    [[nodiscard]]
    std::vector<atom> attached_hydrogens() const;

    [[nodiscard]]
    bool is_a_vdw_candidate() const;

    [[nodiscard]]
    double vdw_radius() const;

    [[nodiscard]]
    double mass() const;

    [[nodiscard]]
    bool is_a_hydrogen() const;

    [[nodiscard]]
    bool is_a_cation() const;

    [[nodiscard]]
    bool is_main_chain() const;

    [[nodiscard]]
    int32_t atom_number() const;

    [[nodiscard]]
    std::string unique_id() const;
};

class ring final : public kdpoint<3>, public aminoacid::component
{
private:
    struct impl;
    std::shared_ptr<impl const> pimpl;

public:
    ring(std::vector<atom> const& atoms, aminoacid const& res);

    ~ring();

    [[nodiscard]]
    std::array<double, 3> const& normal() const;

    [[nodiscard]]
    double radius() const;

    [[nodiscard]]
    bool is_a_pication_candidate() const;

    [[nodiscard]]
    double closest_distance_between_atoms(ring const& other) const;

    [[nodiscard]]
    double angle_between_normals(ring const& other) const;

    [[nodiscard]]
    double angle_between_normal_and_centres_joining(ring const& other) const;

    [[nodiscard]]
    atom atom_closest_to(atom const& atom) const;

    [[nodiscard]]
    std::string name() const;
};

class ionic_group final : public kdpoint<3>, public aminoacid::component
{
private:
    struct impl;
    std::shared_ptr<impl const> pimpl;

public:
    ionic_group(std::vector<atom> const& atoms, int const& charge, aminoacid const& res);

    ~ionic_group();

    [[nodiscard]]
    int charge() const;

    [[nodiscard]]
    double ionion_energy_q() const;

    [[nodiscard]]
    std::string name() const;
};
}