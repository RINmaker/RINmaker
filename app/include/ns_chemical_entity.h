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
    std::shared_ptr<impl> _pimpl;

    aminoacid();

public:
    class component
    {
    private:
        std::weak_ptr<aminoacid::impl> _res_pimpl;

    protected:
        explicit component(aminoacid const& res) : _res_pimpl{res._pimpl}
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

    aminoacid(
        gemmi::Residue const& residue,
        gemmi::Chain const& chain,
        gemmi::Model const& model,
        gemmi::Structure const& protein,
        std::optional<gemmi::Helix> const& helix);

    aminoacid(
        gemmi::Residue const& residue,
        gemmi::Chain const& chain,
        gemmi::Model const& model,
        gemmi::Structure const& protein,
        std::optional<gemmi::Sheet::Strand> const& strand);

    ~aminoacid();

    [[nodiscard]]
    std::vector<atom> const& get_atoms() const;

    [[nodiscard]]
    std::string const& get_protein_name() const;

    [[nodiscard]]
    std::optional<atom> const& get_alpha_carbon() const;

    [[nodiscard]]
    std::optional<atom> const& get_beta_carbon() const;

    [[nodiscard]]
    std::optional<ring> const& get_primary_ring() const;

    [[nodiscard]]
    std::optional<ring> const& get_secondary_ring() const;

    [[nodiscard]]
    std::optional<ionic_group> const& get_positive_ionic_group() const;

    [[nodiscard]]
    std::optional<ionic_group> const& get_negative_ionic_group() const;

    [[nodiscard]]
    std::string const& get_name() const;

    [[nodiscard]]
    std::string const& get_chain_id() const;

    [[nodiscard]]
    std::string const& get_id() const;

    [[nodiscard]]
    int get_sequence_number() const;

    bool operator==(aminoacid const& rhs) const;

    bool operator!=(aminoacid const& rhs) const;

    [[nodiscard]]
    bool satisfies_minimum_sequence_separation(aminoacid const& other, int minimum_separation = cfg::params::seq_sep) const;

    [[nodiscard]]
    explicit operator rin::node() const;

    [[nodiscard]]
    std::string get_secondary_structure_id() const;

    [[nodiscard]]
    std::array<double, 3> const& get_position() const;
};

class atom final : public kdpoint<3>, public aminoacid::component
{
private:
    struct impl;
    std::shared_ptr<impl const> _pimpl;

public:
    atom(gemmi::Atom const& record, aminoacid const& res);

    ~atom();

    [[nodiscard]]
    std::string const& get_name() const;

    [[nodiscard]]
    std::string get_symbol() const;

    [[nodiscard]]
    double get_temp_factor() const;

    [[nodiscard]]
    int get_charge() const;

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
    std::vector<atom> get_attached_hydrogens() const;

    [[nodiscard]]
    bool is_vdw_candidate() const;

    [[nodiscard]]
    double get_vdw_radius() const;

    [[nodiscard]]
    double get_mass() const;

    [[nodiscard]]
    bool is_hydrogen() const;

    [[nodiscard]]
    bool is_cation() const;

    [[nodiscard]]
    bool is_main_chain() const;

    [[nodiscard]]
    int get_atom_number() const;
};

class ring final : public kdpoint<3>, public aminoacid::component
{
private:
    struct impl;
    std::shared_ptr<impl const> _pimpl;

public:
    ring(std::vector<atom> const& atoms, aminoacid const& res);

    ~ring();

    [[nodiscard]]
    std::array<double, 3> const& get_normal() const;

    [[nodiscard]]
    bool is_pication_candidate() const;

    [[nodiscard]]
    double get_distance_between_closest_atoms(ring const& other) const;

    [[nodiscard]]
    double get_angle_between_normals(ring const& other) const;

    [[nodiscard]]
    double get_angle_between_normal_and_centers_joining(ring const& other) const;

    [[nodiscard]]
    std::string get_name() const;
};

class ionic_group final : public kdpoint<3>, public aminoacid::component
{
private:
    struct impl;
    std::shared_ptr<impl const> _pimpl;

public:
    ionic_group(std::vector<atom> const& atoms, int const& charge, aminoacid const& res);

    ~ionic_group();

    [[nodiscard]]
    int get_charge() const;

    [[nodiscard]]
    double get_ionion_energy_q() const;

    [[nodiscard]]
    std::string get_name() const;
};
}