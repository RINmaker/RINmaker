#pragma once

#pragma warning(push, 0)

#include <gemmi/pdb.hpp>

#pragma warning(pop)

#include <string>
#include <vector>
#include <array>
#include <memory>
#include <optional>

#include "config.h"

#include "rin_graph.h"

#include "spatial/kdpoint.h"

namespace chemical_entity
{
class atom;

class ring;

class ionic_group;

class aminoacid
{
private:
    struct impl;
    std::shared_ptr<impl> pimpl;

    aminoacid();

public:
    class component
    {
    private:
        std::weak_ptr<aminoacid::impl> res_impl;

    protected:
        explicit component(aminoacid const& res) : res_impl{res.pimpl}
        {}

    public:
        [[nodiscard]]
        aminoacid get_residue() const;
    };

    friend class aminoacid::component;

public:
    aminoacid(
        gemmi::Residue const& residue,
        gemmi::Chain const& chain,
        gemmi::Model const& model,
        gemmi::Structure const& protein);

    ~aminoacid();

    [[nodiscard]]
    std::vector<atom> const& get_atoms() const;

    [[nodiscard]]
    std::string const& get_protein_name() const;

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

    [[nodiscard]]
    std::string secondary_structure_id() const;

    [[nodiscard]]
    std::array<double, 3> const& position() const;
};

class atom final : public kdpoint<3>, public aminoacid::component
{
private:
    struct impl;
    std::shared_ptr<impl const> pimpl;

public:
    atom(gemmi::Atom const& record, aminoacid const& res);

    ~atom();

    [[nodiscard]]
    std::string const& name() const;

    [[nodiscard]]
    std::string symbol() const;

    [[nodiscard]]
    double temp_factor() const;

    [[nodiscard]]
    int charge() const;

    [[nodiscard]]
    bool in_positive_ionic_group() const;

    [[nodiscard]]
    bool in_negative_ionic_group() const;

    [[nodiscard]]
    bool is_hydrogen_donor() const;

    [[nodiscard]]
    int how_many_hydrogen_can_donate() const;

    [[nodiscard]]
    bool is_hydrogen_acceptor() const;

    [[nodiscard]]
    int how_many_hydrogen_can_accept() const;

    [[nodiscard]]
    std::vector<atom> attached_hydrogens() const;

    [[nodiscard]]
    bool is_vdw_candidate() const;

    [[nodiscard]]
    double vdw_radius() const;

    [[nodiscard]]
    double mass() const;

    [[nodiscard]]
    bool is_hydrogen() const;

    [[nodiscard]]
    bool is_cation() const;

    [[nodiscard]]
    bool is_main_chain() const;

    [[nodiscard]] int32_t atom_number() const;
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