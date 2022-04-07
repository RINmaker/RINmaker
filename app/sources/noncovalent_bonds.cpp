#include "noncovalent_bonds.h"

#include <utility>

#include "chemical_entity.h"

using std::shared_ptr;
using chemical_entity::aminoacid;

bonds::base::base(double length, double energy) : _length(length), _energy(energy)
{}

double bonds::base::get_length() const
{ return _length; }

double bonds::base::get_energy() const
{ return _energy; }

bool bonds::base::operator<(base const& rhs) const
{ return _energy < rhs._energy || (_energy == rhs._energy && _length < rhs._length); }

bool bonds::base::operator>(base const& rhs) const
{ return rhs < *this; }

bonds::computed::computed(aminoacid const& source, aminoacid const& target, double distance, double energy) :
        base(distance, energy),
        _source(source),
        _target(target)
{}

std::string bonds::computed::id() const
{ return prelude::sort(_source.id(), _target.id()); }

chemical_entity::aminoacid const& bonds::computed::source() const
{ return _source; }

chemical_entity::aminoacid const& bonds::computed::target() const
{ return _target; }

bonds::generico::generico(aminoacid const& source, aminoacid const& target)
        : computed(source, target, source.distance(target), 0) // TODO res horribilis
{}

std::string bonds::generico::get_interaction() const
{
    std::string get_interaction = "GENERIC:";
    switch (rin::parameters::global::instance().get().interaction_type())
    {
    case rin::parameters::interaction_type_t::ALPHA_BACKBONE:
        get_interaction += "CA";
        break;

    case rin::parameters::interaction_type_t::BETA_BACKBONE:
        get_interaction += "CB";
        break;

    case rin::parameters::interaction_type_t::NONCOVALENT_BONDS:
        get_interaction += "CLOSEST";
        break;
    }

    return get_interaction;
}

std::string bonds::generico::get_type() const
{ return "generic"; } // TODO config

bonds::generico::operator rin::edge() const
{ return rin::edge(*this); }


//Returns a pair of Sigmaij Epsilonij
std::pair<double, double> bonds::hydrogen::getSigmaEpsilon(
        chemical_entity::atom const& donor, chemical_entity::atom const& acceptor)
{
    std::function<bool(string const&, int, string const&, int)> compare =
            [&](string const& donor_element, int donor_charge, string const& acceptor_element, int acceptor_charge)
            {
                return donor.symbol() == donor_element
                       && donor.charge() == donor_charge
                       && acceptor.symbol() == acceptor_element
                       && acceptor.charge() == acceptor_charge;
            };

    if (compare("N", 0, "N", 0)) return std::make_pair(1.99, -3.00);
    if (compare("N", 0, "O", 0)) return std::make_pair(1.89, -3.50);
    if (compare("O", 0, "N", 0)) return std::make_pair(1.89, -4.00);
    if (compare("O", 0, "O", 0)) return std::make_pair(1.79, -4.25);
    if (compare("N", 1, "N", 0)) return std::make_pair(1.99, -4.50);
    if (compare("N", 1, "O", 0)) return std::make_pair(1.89, -5.25);

    if (compare("N", 0, "O", -1)) return std::make_pair(1.89, -5.25);
    if (compare("N", 1, "O", -1)) return std::make_pair(1.89, -7.00);
    if (compare("O", 0, "O", -1)) return std::make_pair(1.79, -6.375);

    return std::make_pair(0, 0); //TODO log: non dovrebbe accadere
}

double bonds::hydrogen::energy(
        chemical_entity::atom const& donor, chemical_entity::atom const& acceptor, chemical_entity::atom const* hydrogen)
{
    std::pair<double, double> sigmaEpsilon = getSigmaEpsilon(donor, acceptor);
    double sigma = sigmaEpsilon.first;
    double epsilon = sigmaEpsilon.second;
    double distance = hydrogen->distance(acceptor);

    double sigma_distance_12 = pow(sigma / distance, 12);
    double sigma_distance_10 = pow(sigma / distance, 10);

    return 4 * epsilon * (sigma_distance_12 - sigma_distance_10);
}

bonds::hydrogen::hydrogen(
        chemical_entity::atom const& acceptor, chemical_entity::atom const& donor, chemical_entity::atom const* hydrogen, double angle) :
        computed(acceptor.res(), donor.res(), acceptor.distance(donor), energy(donor, acceptor, hydrogen)),
        _acceptor(acceptor),
        _donor(donor),
        _hydrogen(hydrogen),
        _angle(angle)
{}

chemical_entity::atom const& bonds::hydrogen::acceptor() const
{ return _acceptor; }

chemical_entity::atom const& bonds::hydrogen::donor() const
{ return _donor; }

chemical_entity::atom const* bonds::hydrogen::acceptor_ptr() const
{ return &_acceptor; }

chemical_entity::atom const* bonds::hydrogen::donor_ptr() const
{ return &_donor; }

chemical_entity::atom const* bonds::hydrogen::hydrogen_ptr() const
{ return _hydrogen; }

double bonds::hydrogen::get_angle() const
{ return _angle; }

std::string bonds::hydrogen::get_interaction() const
{
    string donorChain = _donor.is_main_chain() ? "MC" : "SC";
    string acceptorChain = _acceptor.is_main_chain() ? "MC" : "SC";

    return "HBOND:" + acceptorChain + "_" + donorChain;
}

bonds::hydrogen::operator rin::edge() const
{ return rin::edge(*this); }

std::string bonds::hydrogen::get_type() const
{ return "hydrogen"; } // TODO va in config

bonds::ionic::ionic(chemical_entity::ionic_group const& negative, chemical_entity::ionic_group const& positive) :
        computed(
                negative.res(), positive.res(), negative.distance(positive), (constant::ion_ion_k * positive.ionion_energy_q() * negative.ionion_energy_q() / (negative.distance(positive)))),
        _negative(negative),
        _positive(positive)
{}

chemical_entity::ionic_group const& bonds::ionic::positive() const
{ return _positive; }

chemical_entity::ionic_group const& bonds::ionic::negative() const
{ return _negative; }

std::string bonds::ionic::get_interaction() const
{ return "IONIC:SC_SC"; }

bonds::ionic::operator rin::edge() const
{ return rin::edge(*this); }

std::string bonds::ionic::get_type() const
{ return "ionic"; }

bonds::pication::pication(chemical_entity::ring const& ring, chemical_entity::atom const& cation, double angle) :
        computed(ring.res(), cation.res(), ring.distance(cation), 9.6),// TODO va in config
        _cation(cation),
        _ring(ring),
        _angle(angle)
{}

chemical_entity::ring const& bonds::pication::ring() const
{ return _ring; }

chemical_entity::atom const& bonds::pication::cation() const
{ return _cation; }

double bonds::pication::angle() const
{ return _angle; }

std::string bonds::pication::get_interaction() const
{ return "PICATION:SC_SC"; }

bonds::pication::operator rin::edge() const
{ return rin::edge(*this); }

std::string bonds::pication::get_type() const
{ return "pication"; }


bonds::pipistack::pipistack(const chemical_entity::ring& source_ring, const chemical_entity::ring& target_ring, double angle) :
        computed(source_ring.res(), target_ring.res(), source_ring.distance(target_ring), 9.6), // TODO va in config
        _source_ring(source_ring),
        _target_ring(target_ring),
        _angle(angle)
{}

chemical_entity::ring const& bonds::pipistack::source_ring() const
{ return _source_ring; }

chemical_entity::ring const& bonds::pipistack::target_ring() const
{ return _target_ring; }

double bonds::pipistack::angle() const
{ return _angle; }

std::string bonds::pipistack::get_interaction() const
{ return "PIPISTACK:SC_SC"; }

bonds::pipistack::operator rin::edge() const
{ return rin::edge(*this); }

std::string bonds::pipistack::get_type() const
{ return "pipistack"; }


bonds::ss::ss(records::ss const& record)
        : base(record.length(), 167), // TODO va in config
          _source_name(record.name_1()),
          _target_name(record.name_2()),
          _source_chain(record.chain_id_1()),
          _target_chain(record.chain_id_2()),
          _source_seq(record.seq_num_1()),
          _target_seq(record.seq_num_2())
{}

std::string bonds::ss::source_id() const
{ return _source_chain + ":" + std::to_string(_source_seq) + ":_:" + _source_name; }

std::string bonds::ss::target_id() const
{ return _target_chain + ":" + std::to_string(_target_seq) + ":_:" + _target_name; }

std::string bonds::ss::get_interaction() const
{ return "SSBOND:SC_SC"; } // TODO config

bonds::ss::operator rin::edge() const
{ return rin::edge(*this); }

std::string bonds::ss::id() const
{ return prelude::sort(source_id(), target_id()); }

std::string bonds::ss::get_type() const
{ return "ss"; } // TODO va in config


double bonds::vdw::energy(chemical_entity::atom const& source_atom, chemical_entity::atom const& target_atom)
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

    // double sigma_distance_6 = pow(sigma / distance, 6);
    double sigma_distance_10 = pow(sigma / distance, 10);

    return 4 * epsilon * (sigma_distance_12 - sigma_distance_10);
}

bonds::vdw::vdw(chemical_entity::atom const& source_atom, chemical_entity::atom const& target_atom) :
        computed(
                source_atom.res(), target_atom.res(), source_atom.distance(target_atom), energy(source_atom, target_atom)),  // TODO config
        _source_atom(source_atom),
        _target_atom(target_atom)
{}

chemical_entity::atom const& bonds::vdw::source_atom() const
{ return _source_atom; }

chemical_entity::atom const& bonds::vdw::target_atom() const
{ return _target_atom; }

std::string bonds::vdw::get_interaction() const
{
    string sourceChain = _source_atom.name() == "C" || _source_atom.name() == "S" ? "MC" : "SC";
    string targetChain = _target_atom.name() == "C" || _target_atom.name() == "S" ? "MC" : "SC";

    return "VDW:" + sourceChain + "_" + targetChain;
}

bonds::vdw::operator rin::edge() const
{ return rin::edge(*this); }

std::string bonds::vdw::get_type() const
{ return "vdw"; }
