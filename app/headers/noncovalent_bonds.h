#pragma once

#include <string>
#include "prelude.h"

#include "rin_graph.h"
#include "pdb_records.h"

#include "energy.h"
#include "rin_params.h"

class network;

namespace chemical_entity
{
class aminoacid;

class atom;

class ring;

class ionic_group;
}

namespace bonds
{
class base
{
private:
    friend class ::network;

    double const _length;
    double const _energy;

protected:
    base(double length, double energy);

public:
    virtual ~base() = default;

public:
    double get_length() const;

    double get_energy() const;

    bool operator<(base const& rhs) const;

    bool operator>(base const& rhs) const;

public:
    virtual std::string get_interaction() const = 0;

    virtual std::string get_type() const = 0;

    [[nodiscard]]
    virtual std::string id() const = 0;

    virtual rin::edge to_edge() const = 0;
};

class computed : public base
{
protected:
    friend class ::network;

    chemical_entity::aminoacid const* _source;
    chemical_entity::aminoacid const* _target;

    computed(chemical_entity::aminoacid const& source, chemical_entity::aminoacid const& target, double distance, double energy);

public:
    virtual ~computed() = default;

    virtual std::string get_type() const = 0;

public:
    chemical_entity::aminoacid const& source() const;

    chemical_entity::aminoacid const& target() const;

    [[nodiscard]]
    std::string id() const override;
};

class generico : public computed
{
private:
    friend class ::network;

    generico(chemical_entity::aminoacid const& source, chemical_entity::aminoacid const& target);

public:
    std::string get_interaction() const;

    std::string get_type() const;

    rin::edge to_edge() const;
};

class hydrogen : public computed
{
private:
    friend class ::network;

    chemical_entity::atom const* _acceptor;
    chemical_entity::atom const* _donor;
    chemical_entity::atom const* _hydrogen;
    double const _angle;

    //Returns a pair of Sigmaij Epsilonij
    std::pair<double, double> getSigmaEpsilon(
            chemical_entity::atom const& donor, chemical_entity::atom const& acceptor);

    double energy(chemical_entity::atom const& donor, chemical_entity::atom const& acceptor, chemical_entity::atom const* hydrogen);

    hydrogen(chemical_entity::atom const& acceptor, chemical_entity::atom const& donor, chemical_entity::atom const* hydrogen, double angle);

public:
    chemical_entity::atom const& acceptor() const;

    chemical_entity::atom const& donor() const;

    chemical_entity::atom const* acceptor_ptr() const;

    chemical_entity::atom const* donor_ptr() const;

    chemical_entity::atom const* hydrogen_ptr() const;

    double get_angle() const;

    std::string get_interaction() const;

    rin::edge to_edge() const;

    std::string get_type() const;
};

class ionic : public computed
{
private:
    friend class ::network;

    chemical_entity::ionic_group const* _negative;
    chemical_entity::ionic_group const* _positive;

    ionic(chemical_entity::ionic_group const& negative, chemical_entity::ionic_group const& positive);

public:
    chemical_entity::ionic_group const& positive() const;

    chemical_entity::ionic_group const& negative() const;

    std::string get_interaction() const;

    rin::edge to_edge() const;

    std::string get_type() const;
};

class pication : public computed
{
private:
    friend class ::network;

    chemical_entity::atom const* _cation;
    chemical_entity::ring const* _ring;
    double _angle;

    pication(chemical_entity::ring const& ring, chemical_entity::atom const& cation, double angle);

public:
    chemical_entity::ring const& ring() const;

    chemical_entity::atom const& cation() const;

    double angle() const;

    std::string get_interaction() const;

    rin::edge to_edge() const;

    std::string get_type() const;
};

class pipistack : public computed
{
private:
    friend class ::network;

    chemical_entity::ring const* _source_ring;
    chemical_entity::ring const* _target_ring;
    double const _angle;

    pipistack(chemical_entity::ring const& source_ring, chemical_entity::ring const& target_ring, double angle);

public:
    chemical_entity::ring const& source_ring() const;

    chemical_entity::ring const& target_ring() const;

    double angle() const;

    std::string get_interaction() const;

    rin::edge to_edge() const;

    std::string get_type() const;
};

class ss : public base
{
private:
    friend class ::network;

    int const _source_seq;
    std::string const _source_name;
    std::string const _source_chain;

    int const _target_seq;
    std::string const _target_chain;
    std::string const _target_name;

    ss(records::ss const& record);

public:
    std::string source_id() const;

    std::string target_id() const;

    std::string get_interaction() const;

    // TODO va in config
    rin::edge to_edge() const;

    std::string id() const;

    std::string get_type() const;
};

class vdw : public computed
{
private:
    friend class ::network;

    chemical_entity::atom const* _source_atom;
    chemical_entity::atom const* _target_atom;

    double energy(chemical_entity::atom const& source_atom, chemical_entity::atom const& target_atom);

    vdw(chemical_entity::atom const& source_atom, chemical_entity::atom const& target_atom);

public:
    chemical_entity::atom const& source_atom() const;

    chemical_entity::atom const& target_atom() const;

    [[nodiscard]]
    std::string get_interaction() const override;

    [[nodiscard]]
    rin::edge to_edge() const override;

    [[nodiscard]]
    std::string get_type() const override;
};
}