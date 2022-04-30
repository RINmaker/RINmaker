#pragma once

#include <string>
#include <memory>

#include "prelude.h"

#include "rin_graph.h"
#include "ns_record.h"

#include "ns_chemical_entity.h"

#include "energy.h"
#include "rin_params.h"

namespace bond
{
class base
{
private:
    double const _length;
    double const _energy;

protected:
    base(double length, double energy);

public:
    virtual ~base() = default;

    [[nodiscard]]
    double get_length() const;

    [[nodiscard]]
    double get_energy() const;

    bool operator<(base const& rhs) const;

    bool operator>(base const& rhs) const;

public:
    [[nodiscard]]
    virtual std::string get_interaction() const = 0;

    [[nodiscard]]
    virtual std::string get_type() const = 0;

    [[nodiscard]]
    virtual std::string get_id() const = 0;

    [[nodiscard]]
    virtual std::string get_id_simple() const = 0;

    [[nodiscard]]
    virtual explicit operator rin::edge() const = 0;
};

class generic_bond final : public base
{
private:
    chemical_entity::atom const _source;
    chemical_entity::atom const _target;

public:
    static std::shared_ptr<generic_bond const> test(rin::parameters const& params, chemical_entity::atom const& a, chemical_entity::atom const& b);

    generic_bond(rin::parameters const& params, chemical_entity::atom const& a, chemical_entity::atom const& b);

    [[nodiscard]]
    chemical_entity::aminoacid source() const;

    [[nodiscard]]
    chemical_entity::aminoacid target() const;

    [[nodiscard]]
    std::string get_interaction() const override;

    [[nodiscard]]
    std::string get_type() const override;

    [[nodiscard]]
    std::string get_id() const override;

    [[nodiscard]]
    std::string get_id_simple() const override;

    [[nodiscard]]
    explicit operator rin::edge() const override
    { return rin::edge(*this); }
};

class hydrogen final : public base
{
private:
    chemical_entity::atom const _acceptor;
    chemical_entity::atom const _donor;
    chemical_entity::atom const _hydrogen;

    double const _angle;

    // returns a pair sigma_ij,epsilon_ij
    static std::pair<double, double> getSigmaEpsilon(chemical_entity::atom const& donor, chemical_entity::atom const& acceptor);

    static double energy(chemical_entity::atom const& donor, chemical_entity::atom const& acceptor, chemical_entity::atom const& hydrogen);

public:
    static std::shared_ptr<hydrogen const> test(rin::parameters const& params, chemical_entity::atom const& acceptor, chemical_entity::atom const& donor);

    hydrogen(chemical_entity::atom const& acceptor, chemical_entity::atom const& donor, chemical_entity::atom const& hydrogen, double angle);

    [[nodiscard]]
    chemical_entity::atom const& acceptor() const;

    [[nodiscard]]
    chemical_entity::atom const& donor() const;

    [[nodiscard]]
    chemical_entity::atom const& source_atom() const
    { return acceptor(); }

    [[nodiscard]]
    chemical_entity::atom const& target_atom() const
    { return donor(); }

    [[nodiscard]]
    chemical_entity::atom const& hydrogen_atom() const;

    [[nodiscard]]
    double get_angle() const;

    [[nodiscard]]
    std::string get_interaction() const override;

    [[nodiscard]]
    std::string get_id() const override;

    [[nodiscard]]
    std::string get_id_simple() const override;

    [[nodiscard]]
    explicit operator rin::edge() const override
    { return rin::edge(*this); }

    [[nodiscard]]
    std::string get_type() const override;
};

class ionic final : public base
{
private:
    chemical_entity::ionic_group const _negative;
    chemical_entity::ionic_group const _positive;

public:
    static std::shared_ptr<ionic const> test(rin::parameters const& params, chemical_entity::ionic_group const& a, chemical_entity::ionic_group const& b);

    ionic(chemical_entity::ionic_group const& negative, chemical_entity::ionic_group const& positive);

    [[nodiscard]]
    chemical_entity::ionic_group const& source_positive() const
    { return _positive; }

    [[nodiscard]]
    chemical_entity::ionic_group const& target_negative() const
    { return _negative; }

    [[nodiscard]]
    std::string get_interaction() const override;

    [[nodiscard]]
    explicit operator rin::edge() const override
    { return rin::edge(*this); }

    [[nodiscard]]
    std::string get_type() const override;

    [[nodiscard]]
    std::string get_id() const override;

    [[nodiscard]]
    std::string get_id_simple() const override;
};

class pication : public base
{
private:
    chemical_entity::atom const _cation;
    chemical_entity::ring const _ring;

    double _angle;

public:
    static std::shared_ptr<pication const> test(rin::parameters const& params, chemical_entity::atom const& cation, chemical_entity::ring const& ring);

    pication(chemical_entity::ring const& ring, chemical_entity::atom const& cation, double angle);

    [[nodiscard]]
    chemical_entity::ring const& source_ring() const
    { return _ring; }

    [[nodiscard]]
    chemical_entity::atom const& target_cation() const
    { return _cation; }

    [[nodiscard]]
    double angle() const;

    [[nodiscard]]
    std::string get_interaction() const override;

    [[nodiscard]]
    explicit operator rin::edge() const override
    { return rin::edge(*this); }

    [[nodiscard]]
    std::string get_type() const override;

    [[nodiscard]]
    std::string get_id() const override;

    [[nodiscard]]
    std::string get_id_simple() const override;
};

class pipistack final : public base
{
private:
    chemical_entity::ring const _source_ring;
    chemical_entity::ring const _target_ring;
    double const _angle;

public:
    static std::shared_ptr<pipistack const> test(rin::parameters const& params, chemical_entity::ring const& a, chemical_entity::ring const& b);

    pipistack(chemical_entity::ring const& a, chemical_entity::ring const& b, double angle);

    static double energy(double angle);

    [[nodiscard]]
    chemical_entity::ring const& source_ring() const
    { return _source_ring; }

    [[nodiscard]]
    chemical_entity::ring const& target_ring() const
    { return _target_ring; }

    [[nodiscard]]
    double angle() const;

    [[nodiscard]]
    std::string get_interaction() const override;

    [[nodiscard]]
    explicit operator rin::edge() const override
    { return rin::edge(*this); }

    [[nodiscard]]
    std::string get_type() const override;

    [[nodiscard]]
    std::string get_id() const override;

    [[nodiscard]]
    std::string get_id_simple() const override;
};

class ss final : public base
{
private:
    int const _source_seq;
    std::string const _source_name;
    std::string const _source_chain;

    int const _target_seq;
    std::string const _target_chain;
    std::string const _target_name;

public:
    explicit ss(record::ss const& record);

    [[nodiscard]]
    std::string source_id() const;

    [[nodiscard]]
    std::string target_id() const;

    [[nodiscard]]
    std::string get_interaction() const override;

    [[nodiscard]]
    explicit operator rin::edge() const override
    { return rin::edge(*this); }

    [[nodiscard]]
    std::string get_id() const override;

    [[nodiscard]]
    std::string get_id_simple() const override;

    [[nodiscard]]
    std::string get_type() const override;
};

class vdw final : public base
{
private:
    chemical_entity::atom const _source_atom;
    chemical_entity::atom const _target_atom;

    static double energy(chemical_entity::atom const& source_atom, chemical_entity::atom const& target_atom);

public:
    static std::shared_ptr<vdw const> test(rin::parameters const& params, chemical_entity::atom const& a, chemical_entity::atom const& b);

    vdw(chemical_entity::atom const& a, chemical_entity::atom const& b);

    [[nodiscard]]
    chemical_entity::atom const& source_atom() const
    { return _source_atom; }

    [[nodiscard]]
    chemical_entity::atom const& target_atom() const
    { return _target_atom; }

    [[nodiscard]]
    std::string get_interaction() const override;

    [[nodiscard]]
    explicit operator rin::edge() const override
    { return rin::edge(*this); }

    [[nodiscard]]
    std::string get_type() const override;

    [[nodiscard]]
    std::string get_id() const override;

    [[nodiscard]]
    std::string get_id_simple() const override;
};
}