#pragma once

#include <string>
#include "rin.h"
#include "utils/prelude.h"
#include "entity.h"
#include "parameters.h"
#include "energy.h"

class network;

// legami non covalenti tra due amminoacidi
//
namespace bonds
{
// base class
//
class base
{
private:
    friend class ::network;

    double const _length;
    double const _energy;

protected:
    base(double length, double energy) : _length(length), _energy(energy)
    {}

public:
    virtual ~base() = default;

public:
    double get_length() const
    { return _length; }

    double get_energy() const
    { return _energy; }

    // TODO: sprintf di prelude � stranaprelude::sprintf("[type:%s][id:%s][interaction:%s]", get_type(), get_key(), get_interaction())
    // std::string pretty() const { return prelude::sprintf("[type:%s][id:%s][interaction:%s]", get_type(), get_key(), get_interaction()); }

    // confronta prima l'energia "per tipo di legame"; se uguale, va per distanza.
    //
    bool operator<(base const& rhs) const
    { return _energy < rhs._energy || (_energy == rhs._energy && _length < rhs._length); }

    bool operator>(base const& rhs) const
    { return rhs < *this; }

public:
    virtual std::string get_interaction() const = 0;

    virtual std::string get_type() const = 0;

    [[nodiscard]]
    virtual std::string id() const = 0;

    virtual rin::edge to_edge() const = 0;
};

// legami che vengono computati di fatto (esempio negativo: l'ss bond � solo parsato dal pdb)
//
class computed : public base
{
protected:
    friend class ::network;

    entities::aminoacid const* _source;
    entities::aminoacid const* _target;

    computed(entities::aminoacid const& source, entities::aminoacid const& target, double distance, double energy)
            : base(distance, energy), _source(&source), _target(&target)
    {}

public:
    virtual ~computed() = default;

    virtual std::string get_type() const = 0;

public:
    entities::aminoacid const& source() const
    { return *_source; }

    entities::aminoacid const& target() const
    { return *_target; }

    [[nodiscard]]
    std::string id() const override
    { return prelude::sort(_source->id(), _target->id()); }
};

// legame generico tra carbonio alpha o beta
//
class generico : public computed
{
private:
    friend class ::network;

    generico(entities::aminoacid const& source, entities::aminoacid const& target)
            : computed(source, target, source.distance(target), 0) // TODO res horribilis
    {}

public:
    std::string get_interaction() const
    {
        std::string get_interaction = "GENERIC:";
        switch (parameters::get_net_policy())
        {
            case parameters::policy::CA:
                get_interaction += "CA";
                break;

            case parameters::policy::CB:
                get_interaction += "CB";
                break;

            case parameters::policy::CLOSEST:
                get_interaction += "CLOSEST";
                break;
        }

        return get_interaction;
    }

    std::string get_type() const
    { return "generic"; } // TODO va in config
    rin::edge to_edge() const
    { return rin::edge(*this); }
};

// legame idrogeno
//
class hydrogen : public computed
{
private:
    friend class ::network;

    entities::atom const* _acceptor;
    entities::atom const* _donor;
    entities::atom const* _hydrogen;
    double const _angle;

    //Returns a pair of Sigmaij Epsilonij
    std::pair<double, double> getSigmaEpsilon(entities::atom const& donor,
                                          entities::atom const& acceptor)
    {
        std::function<bool(string const&, int, string const&, int)> compare =
          [&](string const& donor_element, int donor_charge, string const& acceptor_element, int acceptor_charge)
        {
            return donor.symbol() == donor_element
                   && donor.charge() == donor_charge
                   && acceptor.symbol() == acceptor_element
                   && acceptor.charge() == acceptor_charge;
        };

        if (compare("N", 0, "N", 0)) return  std::make_pair(1.99, -3.00);
        if (compare("N", 0, "O", 0)) return  std::make_pair(1.89, -3.50);
        if (compare("O", 0, "N", 0)) return  std::make_pair(1.89, -4.00);
        if (compare("O", 0, "O", 0)) return  std::make_pair(1.79, -4.25);
        if (compare("N", 1, "N", 0)) return  std::make_pair(1.99, -4.50);
        if (compare("N", 1, "O", 0)) return  std::make_pair(1.89, -5.25);

        if (compare("N", 0, "O", -1)) return  std::make_pair(1.89, -5.25);
        if (compare("N", 1, "O", -1)) return  std::make_pair(1.89, -7.00);
        if (compare("O", 0, "O", -1)) return  std::make_pair(1.79, -6.375);

        return std::make_pair(0, 0); //TODO log: non dovrebbe accadere
    }

    double energy(entities::atom const& donor, entities::atom const& acceptor, entities::atom const* hydrogen)
    {
        std::pair<double, double> sigmaEpsilon = getSigmaEpsilon(donor, acceptor);
        double sigma = sigmaEpsilon.first;
        double epsilon = sigmaEpsilon.second;
        double distance = hydrogen->distance(acceptor);

        double sigma_distance_12 = pow(sigma / distance, 12);
        double sigma_distance_10 = pow(sigma / distance, 10);

        return 4 * epsilon * (sigma_distance_12 - sigma_distance_10);
    }

    hydrogen(entities::atom const& acceptor, entities::atom const& donor, entities::atom const* hydrogen, double angle)
            : computed(acceptor.res(), donor.res(), acceptor.distance(donor), energy(donor, acceptor, hydrogen)),
            _acceptor(&acceptor), _donor(&donor), _hydrogen(hydrogen), _angle(angle)
    {}

public:
    entities::atom const& acceptor() const
    { return *_acceptor; }
    entities::atom const& donor() const
    { return *_donor; }

    entities::atom const* acceptor_ptr() const
    { return _acceptor; }
    entities::atom const* donor_ptr() const
    { return _donor; }
    entities::atom const* hydrogen_ptr() const
    { return _hydrogen; }

    double get_angle() const
    { return _angle; }

    std::string get_interaction() const
    {
        entities::atom const& donor = *_donor;
        entities::atom const& acceptor = *_acceptor;

        string donorChain    = donor.is_main_chain() ? "MC" : "SC";
        string acceptorChain = acceptor.is_main_chain() ? "MC" : "SC";

        return "HBOND:" + acceptorChain + "_" + donorChain;
    }

    rin::edge to_edge() const
    { return rin::edge(*this); }

    std::string get_type() const
    { return "hydrogen"; } // TODO va in config
};

// legame tra gruppi ionici
//
class ionic : public computed
{
private:
    friend class ::network;

    entities::ionic_group const* _negative;
    entities::ionic_group const* _positive;

    ionic(entities::ionic_group const& negative, entities::ionic_group const& positive)
            : computed(
            negative.res(),
            positive.res(),
            negative.distance(positive),
            (constant::ion_ion_k * positive.ionion_energy_q() * negative.ionion_energy_q() / (negative.distance(positive)))
    )
            , _negative(&negative), _positive(&positive)
    {}

public:
    entities::ionic_group const& positive() const
    { return *_positive; }

    entities::ionic_group const& negative() const
    { return *_negative; }

    std::string get_interaction() const
    { return "IONIC:SC_SC"; }

    rin::edge to_edge() const
    { return rin::edge(*this); }

    std::string get_type() const
    { return "ionic"; }
};

// legame pi-catione
//
class pication : public computed
{
private:
    friend class ::network;

    entities::atom const* _cation;
    entities::ring const* _ring;
    double _angle;

    pication(entities::ring const& ring, entities::atom const& cation, double angle)
            : computed(ring.res(), cation.res(), ring.distance(cation), 9.6) // TODO va in config
            , _cation(&cation), _ring(&ring), _angle(angle)
    {}

public:
    entities::ring const& ring() const
    { return *_ring; }

    entities::atom const& cation() const
    { return *_cation; }

    double angle() const
    { return _angle; }

    std::string get_interaction() const
    { return "PICATION:SC_SC"; }
    rin::edge to_edge() const
    { return rin::edge(*this); }

    std::string get_type() const
    { return "pication"; }
};

// legame pipistack
//
class pipistack : public computed
{
private:
    friend class ::network;

    entities::ring const* _source_ring;
    entities::ring const* _target_ring;
    double const _angle;

    pipistack(entities::ring const& source_ring, entities::ring const& target_ring, double angle)
            : computed(source_ring.res(), target_ring.res(), source_ring.distance(target_ring), 9.6) // TODO va in config
            , _source_ring(&source_ring), _target_ring(&target_ring), _angle(angle)
    {}

public:
    entities::ring const& source_ring() const
    { return *_source_ring; }

    entities::ring const& target_ring() const
    { return *_target_ring; }

    double angle() const
    { return _angle; }

    std::string get_interaction() const
    { return "PIPISTACK:SC_SC"; }
    rin::edge to_edge() const
    { return rin::edge(*this); }

    std::string get_type() const
    { return "pipistack"; }
};

// ponte ss
//
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

    ss(records::ss const& record)
            : base(record.length(), 167) // TODO va in config

            ,
            _source_name(record.name_1()),
            _target_name(record.name_2()),
            _source_chain(record.chain_id_1()),
            _target_chain(record.chain_id_2()),
            _source_seq(record.seq_num_1()),
            _target_seq(record.seq_num_2())
    {}

public:
    std::string source_id() const
    { return _source_chain + ":" + std::to_string(_source_seq) + ":_:" + _source_name; }

    std::string target_id() const
    { return _target_chain + ":" + std::to_string(_target_seq) + ":_:" + _target_name; }

    std::string get_interaction() const
    { return "SSBOND:SC_SC"; } // TODO va in config
    rin::edge to_edge() const
    { return rin::edge(*this); }

    std::string id() const
    { return prelude::sort(source_id(), target_id()); }

    std::string get_type() const
    { return "ss"; } // TODO va in config
};

// legame di van der waals
//
class vdw : public computed
{
private:
    friend class ::network;

    entities::atom const* _source_atom;
    entities::atom const* _target_atom;

    double energy(entities::atom const& source_atom, entities::atom const& target_atom)
    {
        double* source_opts = get_vdw_opsl_values(source_atom.res().name(), source_atom.name(), source_atom.symbol());
        double* target_opts = get_vdw_opsl_values(target_atom.res().name(), target_atom.name(), target_atom.symbol());

        double source_sigma = source_opts[1];
        double target_sigma = target_opts[1];
        double source_epsilon = source_opts[2];
        double target_epsilon = target_opts[2];

        double sigma = sqrt(source_sigma * target_sigma);
        double epsilon = sqrt(source_epsilon * target_epsilon);
        double distance = source_atom.distance(target_atom); //Rij is the distance between center

        double sigma_distance_12 = pow(sigma / distance, 12);
        double sigma_distance_6 = pow(sigma / distance, 6);

        return 4 * epsilon * (sigma_distance_12 - sigma_distance_6);
    }

    vdw(entities::atom const& source_atom, entities::atom const& target_atom)
            : computed(source_atom.res(), target_atom.res(), source_atom.distance(target_atom), energy(source_atom, target_atom))  // TODO va in config
            , _source_atom(&source_atom), _target_atom(&target_atom)
    {}

public:
    entities::atom const& source_atom() const
    { return *_source_atom; }

    entities::atom const& target_atom() const
    { return *_target_atom; }

    std::string get_interaction() const
    {
        entities::atom const& source = *_source_atom;
        entities::atom const& target = *_target_atom;

        string sourceChain = source.name() == "C" || source.name() == "S" ? "MC" : "SC";
        string targetChain = target.name() == "C" || target.name() == "S" ? "MC" : "SC";

        return "VDW:" + sourceChain + "_" + targetChain;
    }

    rin::edge to_edge() const
    { return rin::edge(*this); }

    std::string get_type() const
    { return "vdw"; }
};
}