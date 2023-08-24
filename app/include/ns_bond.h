#pragma once

#pragma warning(push, 0)

#include <cstdint> // apparently, gemmi does not include it properly (but it needs it).
#include <gemmi/pdb.hpp>

#pragma warning(pop)

#include <string>
#include <memory>

#include "rin_graph.h"

#include "ns_chemical_entity.h"

#include "energy.h"
#include "rin_params.h"

namespace bond
{
class base
{
protected:
    double _length;
    double _energy;

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
    virtual std::string get_id() const = 0;

    [[nodiscard]]
    virtual std::string get_id_simple() const = 0;

    // returns the id of the pair of aminoacids
    [[nodiscard]]
    virtual std::string get_pair_id() const = 0;

    [[nodiscard]]
    virtual explicit operator rin::edge() const = 0;
};

class generic_bond : public base
{
private:
    chemical_entity::atom const _source;
    chemical_entity::atom const _target;

protected:
    generic_bond(chemical_entity::atom const& a, chemical_entity::atom const& b, double distance);

public:
    [[nodiscard]]
    chemical_entity::aminoacid get_source() const;

    [[nodiscard]]
    chemical_entity::aminoacid get_target() const;

    [[nodiscard]]
    std::string get_interaction() const override;

    [[nodiscard]]
    std::string get_id() const override;

    [[nodiscard]]
    std::string get_id_simple() const override;

    [[nodiscard]]
    std::string get_pair_id() const override;

    [[nodiscard]]
    explicit operator rin::edge() const override
    { return rin::edge(*this); }

    [[nodiscard]]
    virtual bool has_energy() const = 0;
};

class hydrophobic final : public generic_bond
{
public:
    static std::shared_ptr<hydrophobic const> test(rin::parameters const& params, chemical_entity::atom const& a, chemical_entity::atom const& b);

    hydrophobic(chemical_entity::atom const& c1, chemical_entity::atom const& c2);

    [[nodiscard]]
    std::string get_interaction() const override;

    [[nodiscard]]
    bool has_energy() const override
    { return true; }
};

class contact final : public generic_bond
{
public:
    static std::shared_ptr<contact const> test(rin::parameters const& params, chemical_entity::atom const& a, chemical_entity::atom const& b);

    contact(chemical_entity::atom const& a, chemical_entity::atom const& b);

    [[nodiscard]]
    bool has_energy() const override
    { return false; }
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
    chemical_entity::atom const& get_acceptor() const;

    [[nodiscard]]
    chemical_entity::atom const& get_donor() const;

    [[nodiscard]]
    chemical_entity::atom const& get_source_atom() const
    { return get_acceptor(); }

    [[nodiscard]]
    chemical_entity::atom const& get_target_atom() const
    { return get_donor(); }

    [[nodiscard]]
    chemical_entity::atom const& get_hydrogen_atom() const;

    [[nodiscard]]
    double get_angle() const;

    [[nodiscard]]
    std::string get_interaction() const override;

    [[nodiscard]]
    std::string get_id() const override;

    [[nodiscard]]
    std::string get_id_simple() const override;

    [[nodiscard]]
    std::string get_pair_id() const override;

    [[nodiscard]]
    explicit operator rin::edge() const override
    { return rin::edge(*this); }
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
    chemical_entity::ionic_group const& get_source_positive() const
    { return _positive; }

    [[nodiscard]]
    chemical_entity::ionic_group const& get_target_negative() const
    { return _negative; }

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
    std::string get_pair_id() const override;
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
    chemical_entity::ring const& get_source_ring() const
    { return _ring; }

    [[nodiscard]]
    chemical_entity::atom const& get_target_cation() const
    { return _cation; }

    [[nodiscard]]
    double get_angle() const;

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
    std::string get_pair_id() const override;

    static double getKappa(const chemical_entity::atom &cation);

    static double getAlpha(const chemical_entity::ring &ring);

    static double energy(const chemical_entity::ring &ring, const chemical_entity::atom &cation);
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
    chemical_entity::ring const& get_source_ring() const
    { return _source_ring; }

    [[nodiscard]]
    chemical_entity::ring const& get_target_ring() const
    { return _target_ring; }

    [[nodiscard]]
    double get_angle() const;

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
    std::string get_pair_id() const override;
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
    explicit ss(gemmi::Connection const& connection);

    [[nodiscard]]
    std::string get_source_id() const;

    [[nodiscard]]
    std::string get_target_id() const;

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
    std::string get_pair_id() const override;
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
    chemical_entity::atom const& get_source_atom() const
    { return _source_atom; }

    [[nodiscard]]
    chemical_entity::atom const& get_target_atom() const
    { return _target_atom; }

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
    std::string get_pair_id() const override;
};
}