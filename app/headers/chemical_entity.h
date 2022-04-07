#pragma once

#include <string>
#include <vector>
#include <array>
#include <tuple>

#include "config.h"
#include "prelude.h"

#include "log_manager.h"

#include "rin_graph.h"
#include "pdb_records.h"

#include "spatial/geometry.h"
#include "spatial/kdpoint.h"

namespace rin
{
struct maker;
}

namespace structure
{
class base;
}

namespace chemical_entity
{

class atom;

class ring;

class ionic_group;

class aminoacid : public kdpoint<3>
{
private:
    std::vector<atom const*> _atoms = std::vector<atom const*>();
    atom const* _alpha_carbon = nullptr;
    atom const* _beta_carbon = nullptr;

    ring const* _primary_ring = nullptr;
    ring const* _secondary_ring = nullptr;

    ionic_group const* _positive_ionic_group = nullptr;
    ionic_group const* _negative_ionic_group = nullptr;

private:
    std::string _chain_id;

    std::string _name;

    int _sequence_number = 0;

    std::string _id;
    structure::base* _secondary_structure;

private:
    friend struct rin::maker;

    std::string _pdb_name;

    explicit aminoacid(std::vector<records::atom> const& records, std::string pdb_name);

    ~aminoacid();

public:
    [[nodiscard]]
    std::vector<atom const*> const& atoms() const
    { return _atoms; }

    [[nodiscard]]
    std::string const& pdb_name() const
    { return _pdb_name; }

    [[nodiscard]]
    atom const* ca() const
    { return _alpha_carbon; }

    [[nodiscard]]
    atom const* cb() const
    { return _beta_carbon; }

    [[nodiscard]]
    ring const* primary_ring() const
    { return _primary_ring; }

    [[nodiscard]]
    ring const* secondary_ring() const
    { return _secondary_ring; }

    [[nodiscard]]
    ionic_group const* positive_ionic_group() const
    { return _positive_ionic_group; }

    [[nodiscard]]
    ionic_group const* negative_ionic_group() const
    { return _negative_ionic_group; }

public:
    [[nodiscard]]
    std::string const& name() const
    { return _name; }

    [[nodiscard]]
    std::string const& chain_id() const
    { return _chain_id; }

    [[nodiscard]]
    std::string const& id() const
    { return _id; }

    [[nodiscard]]
    int sequence_number() const
    { return _sequence_number; }

public:
    bool operator==(aminoacid const& rhs) const
    { return _id == rhs._id; }

    bool operator!=(aminoacid const& rhs) const
    { return !(*this == rhs); }

    [[nodiscard]]
    bool
    satisfies_minimum_separation(aminoacid const& aa, int minimum_separation = cfg::params::seq_sep) const
    {
        if (*this == aa)
        {
            return false;
        }

        if (_chain_id != aa._chain_id)
        {
            return true;
        }

        return abs(_sequence_number - aa._sequence_number) >= minimum_separation;
    }

    [[nodiscard]]
    rin::node to_node() const
    { return rin::node(*this); } // TODO può diventare un operatore di conversione

public:
    void make_secondary_structure();

    void make_secondary_structure(records::helix const& record);

    void make_secondary_structure(records::sheet_piece const& record);

    [[nodiscard]]
    std::string secondary_structure_id() const;
};

class component
{
protected:
    friend class aminoacid;

    aminoacid const& _res;
    explicit component(aminoacid const& res) : _res(res)
    {}

    virtual ~component() = default;

public:
    [[nodiscard]]
    aminoacid const& res() const
    { return _res; }
};

class atom final : public kdpoint<3>, public component
{
private:
    friend class aminoacid;

    atom(records::atom const& record, aminoacid const& res)
            : kdpoint<3>({record.x(), record.y(), record.z()}), component(res), _record(record)
    {}

    ~atom() override = default;

private:
    // the record it was parsed from
    records::atom _record;

public:
    [[nodiscard]]
    std::string const& name() const
    { return _record.name(); }

    [[nodiscard]]
    std::string const& symbol() const
    { return _record.element_name(); }

    [[nodiscard]]
    double temp_factor() const
    { return _record.temp_factor(); }

    [[nodiscard]]
    int charge() const
    { return _record.charge(); }

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
    bool is_a_hydrogen() const
    { return symbol() == "H"; }

    [[nodiscard]]
    bool is_a_cation() const;

    [[nodiscard]]
    bool is_main_chain() const
    {
        return
                this->name() == "C" ||
                this->name() == "O" ||
                this->name() == "H" ||
                this->name() == "HA" ||
                this->name() == "N";
    }
};

class ring final : public kdpoint<3>, public component
{
private:
    std::vector<atom const*> _atoms;

    std::array<double, 3> _normal{};
    double _mean_radius;

private:
    friend class aminoacid;

    ring(std::vector<atom const*> const& atoms, aminoacid const& res);

    ~ring() override = default;

public:
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

class ionic_group final : public kdpoint<3>, public component
{
private:
    std::vector<atom const*> const _atoms;
    int const _charge;

private:
    friend class aminoacid;

    ionic_group(std::vector<atom const*> const& atoms, int const& charge, aminoacid const& res);

    ~ionic_group() override = default;

public:
    [[nodiscard]]
    int charge() const
    { return _charge; }

    [[nodiscard]]
    double ionion_energy_q() const;

    [[nodiscard]]
    string name() const;
};
}

string getNameFromAtoms(std::vector<const chemical_entity::atom*> const& atoms, string const& delimiter = ":");
