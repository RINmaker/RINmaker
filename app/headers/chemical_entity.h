#pragma once

#include <string>
#include <vector>
#include <array>
#include <tuple>
#include <memory>

#include "config.h"
#include "prelude.h"

#include "log_manager.h"

#include "rin_graph.h"
#include "pdb_records.h"

#include "secondary_structures.h"

#include "spatial/geometry.h"
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
    impl* pimpl;

public:
    class component
    {
    protected:
        aminoacid const& _res;

        explicit component(aminoacid const& res) : _res(res)
        {}

    public:
        [[nodiscard]]
        aminoacid const& res() const
        { return _res; }
    };

    aminoacid(std::vector<records::atom> const& records, std::string pdb_name);

    ~aminoacid();

    [[nodiscard]]
    std::vector<atom const*> atoms() const;

    [[nodiscard]]
    std::string const& pdb_name() const;

    [[nodiscard]]
    atom const* ca() const;

    [[nodiscard]]
    atom const* cb() const;

    [[nodiscard]]
    ring const* primary_ring() const;

    [[nodiscard]]
    ring const* secondary_ring() const;

    [[nodiscard]]
    ionic_group const* positive_ionic_group() const;

    [[nodiscard]]
    ionic_group const* negative_ionic_group() const;

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

    void make_secondary_structure();

    void make_secondary_structure(records::helix const& record);

    void make_secondary_structure(records::sheet_piece const& record);

    [[nodiscard]]
    std::string secondary_structure_id() const;
};

class atom final : public kdpoint<3>, public aminoacid::component
{
private:
    struct impl;
    impl* pimpl;

public:
    atom(records::atom const& record, aminoacid const& res);

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
    std::vector<atom const*> attached_hydrogens() const;

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
};

class ring final : public kdpoint<3>, public aminoacid::component
{
private:
    std::vector<atom const*> _atoms;

    std::array<double, 3> _normal{};
    double _mean_radius;

public:
    ring(std::vector<atom const*> const& atoms, aminoacid const& res);

    [[nodiscard]]
    std::array<double, 3> const& normal() const
    { return _normal; }

    [[nodiscard]]
    double radius() const
    { return _mean_radius; }

    [[nodiscard]]
    bool is_a_pication_candidate() const
    {
        string name = res().name();
        return name == "PHE" || name == "TYR" || (name == "TRP" && _atoms.size() == 6);
    }

    [[nodiscard]]
    double closest_distance_between_atoms(ring const& other) const
    {
        double minimum = _atoms[0]->distance(*other._atoms[0]);
        for (auto* atom_1: _atoms)
        {
            for (auto* atom_2: other._atoms)
            {
                double current = atom_1->distance(*atom_2);
                if (current < minimum)
                {
                    minimum = current;
                }
            }
        }

        return minimum;
    }

    [[nodiscard]]
    double angle_between_normals(ring const& other) const
    { return geom::d_angle<3>(_normal, other._normal); }

    [[nodiscard]]
    double angle_between_normal_and_centres_joining(ring const& other) const
    {
        std::array<double, 3> const centres_joining((std::array<double, 3>) (*this - other));
        return geom::d_angle<3>(_normal, centres_joining);
    }

    [[nodiscard]]
    atom const& atom_closest_to(atom const& atom) const;

    [[nodiscard]]
    string name() const;
};

class ionic_group final : public kdpoint<3>, public aminoacid::component
{
private:
    struct impl;
    impl* pimpl;

public:
    ionic_group(std::vector<atom const*> const& atoms, int const& charge, aminoacid const& res);

    ~ionic_group();

    [[nodiscard]]
    int charge() const;

    [[nodiscard]]
    double ionion_energy_q() const;

    [[nodiscard]]
    string name() const;
};
}

string getNameFromAtoms(std::vector<const chemical_entity::atom*> const& atoms, string const& delimiter = ":");
