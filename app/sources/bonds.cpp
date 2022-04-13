#include "bonds.h"

#include <utility>

#include "chemical_entity.h"
#include "rin_network.h"

using std::shared_ptr;
using chemical_entity::aminoacid;

bond::base::base(double length, double energy) : _length(length), _energy(energy)
{}

double bond::base::get_length() const
{ return _length; }

double bond::base::get_energy() const
{ return _energy; }

bool bond::base::operator<(base const& rhs) const
{ return _energy < rhs._energy || (_energy == rhs._energy && _length < rhs._length); }

bool bond::base::operator>(base const& rhs) const
{ return rhs < *this; }

bond::computed::computed(rin::parameters const& params, aminoacid const& source, aminoacid const& target, double distance, double energy) :
        base(distance, energy),
        _params(params),
        _source(source),
        _target(target)
{}

std::string bond::computed::id() const
{ return prelude::concat_lexicographically(_source.id(), _target.id()); }

chemical_entity::aminoacid const& bond::computed::source() const
{ return _source; }

chemical_entity::aminoacid const& bond::computed::target() const
{ return _target; }

bond::generico::generico(rin::parameters const& params, aminoacid const& source, aminoacid const& target)
        : computed(params, source, target, source.distance(target), 0) // TODO res horribilis
{}

std::string bond::generico::get_interaction() const
{
    std::string get_interaction = "GENERIC:";
    switch (_params.interaction_type())
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

std::string bond::generico::get_type() const
{ return "generic"; } // TODO config

//Returns a pair of Sigmaij Epsilonij
std::pair<double, double> bond::hydrogen::getSigmaEpsilon(
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

    return std::make_pair(0, 0); //TODO log: it should not happen
}

double bond::hydrogen::energy(
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

bond::hydrogen::hydrogen(
        rin::parameters const& params, chemical_entity::atom const& acceptor, chemical_entity::atom const& donor, chemical_entity::atom const* hydrogen, double angle) :
        computed(params, acceptor.res(), donor.res(), acceptor.distance(donor), energy(donor, acceptor, hydrogen)),
        _acceptor(acceptor),
        _donor(donor),
        _hydrogen(hydrogen),
        _angle(angle)
{}

chemical_entity::atom const& bond::hydrogen::acceptor() const
{ return _acceptor; }

chemical_entity::atom const& bond::hydrogen::donor() const
{ return _donor; }

chemical_entity::atom const* bond::hydrogen::acceptor_ptr() const
{ return &_acceptor; }

chemical_entity::atom const* bond::hydrogen::donor_ptr() const
{ return &_donor; }

chemical_entity::atom const* bond::hydrogen::hydrogen_ptr() const
{ return _hydrogen; }

double bond::hydrogen::get_angle() const
{ return _angle; }

std::string bond::hydrogen::get_interaction() const
{
    string donorChain = _donor.is_main_chain() ? "MC" : "SC";
    string acceptorChain = _acceptor.is_main_chain() ? "MC" : "SC";

    return "HBOND:" + acceptorChain + "_" + donorChain;
}

std::string bond::hydrogen::get_type() const
{ return "hydrogen"; } // TODO config

bond::ionic::ionic(rin::parameters const& params, chemical_entity::ionic_group const& negative, chemical_entity::ionic_group const& positive) :
        computed(
                params, negative.res(), positive.res(), negative.distance(positive), (constant::ion_ion_k * positive.ionion_energy_q() * negative.ionion_energy_q() / (negative.distance(positive)))),
        _negative(negative),
        _positive(positive)
{}

chemical_entity::ionic_group const& bond::ionic::positive() const
{ return _positive; }

chemical_entity::ionic_group const& bond::ionic::negative() const
{ return _negative; }

std::string bond::ionic::get_interaction() const
{ return "IONIC:SC_SC"; }

std::string bond::ionic::get_type() const
{ return "ionic"; }

bond::pication::pication(rin::parameters const& params, chemical_entity::ring const& ring, chemical_entity::atom const& cation, double angle) :
        computed(params, ring.res(), cation.res(), ring.distance(cation), 9.6),// TODO add formula
        _cation(cation),
        _ring(ring),
        _angle(angle)
{}

chemical_entity::ring const& bond::pication::ring() const
{ return _ring; }

chemical_entity::atom const& bond::pication::cation() const
{ return _cation; }

double bond::pication::angle() const
{ return _angle; }

std::string bond::pication::get_interaction() const
{ return "PICATION:SC_SC"; }

std::string bond::pication::get_type() const
{ return "pication"; }


bond::pipistack::pipistack(rin::parameters const& params, chemical_entity::ring const& source_ring, chemical_entity::ring const& target_ring, double angle) :
        computed(params, source_ring.res(), target_ring.res(), source_ring.distance(target_ring), 9.6), // TODO add formula
        _source_ring(source_ring),
        _target_ring(target_ring),
        _angle(angle)
{}

chemical_entity::ring const& bond::pipistack::source_ring() const
{ return _source_ring; }

chemical_entity::ring const& bond::pipistack::target_ring() const
{ return _target_ring; }

double bond::pipistack::angle() const
{ return _angle; }

std::string bond::pipistack::get_interaction() const
{ return "PIPISTACK:SC_SC"; }

std::string bond::pipistack::get_type() const
{ return "pipistack"; }

bond::ss::ss(records::ss const& record)
        : base(record.length(), 167), // TODO config
          _source_name(record.name_1()),
          _target_name(record.name_2()),
          _source_chain(record.chain_id_1()),
          _target_chain(record.chain_id_2()),
          _source_seq(record.seq_num_1()),
          _target_seq(record.seq_num_2())
{}

std::string bond::ss::source_id() const
{ return _source_chain + ":" + std::to_string(_source_seq) + ":_:" + _source_name; }

std::string bond::ss::target_id() const
{ return _target_chain + ":" + std::to_string(_target_seq) + ":_:" + _target_name; }

std::string bond::ss::get_interaction() const
{ return "SSBOND:SC_SC"; } // TODO config

std::string bond::ss::id() const
{ return prelude::concat_lexicographically(source_id(), target_id()); }

std::string bond::ss::get_type() const
{ return "ss"; }


double bond::vdw::energy(chemical_entity::atom const& source_atom, chemical_entity::atom const& target_atom)
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

bond::vdw::vdw(rin::parameters const& params, chemical_entity::atom const& source_atom, chemical_entity::atom const& target_atom) :
        computed(
                params, source_atom.res(), target_atom.res(), source_atom.distance(target_atom), energy(source_atom, target_atom)),  // TODO config
        _source_atom(source_atom),
        _target_atom(target_atom)
{}

chemical_entity::atom const& bond::vdw::source_atom() const
{ return _source_atom; }

chemical_entity::atom const& bond::vdw::target_atom() const
{ return _target_atom; }

std::string bond::vdw::get_interaction() const
{
    string sourceChain = _source_atom.name() == "C" || _source_atom.name() == "S" ? "MC" : "SC";
    string targetChain = _target_atom.name() == "C" || _target_atom.name() == "S" ? "MC" : "SC";

    return "VDW:" + sourceChain + "_" + targetChain;
}

std::string bond::vdw::get_type() const
{ return "vdw"; }


bool bond::hydrogen::test(network& net, rin::parameters const& params, chemical_entity::atom const& acceptor, chemical_entity::atom const& donor)
{
    if (acceptor.res().satisfies_minimum_separation(donor.res()))
    {
        if (!(acceptor.res() == donor.res()))
        {
            auto hydrogens = donor.attached_hydrogens();
            for (auto* h: hydrogens)
            {
                std::array<double, 3> const da = (std::array<double, 3>) (acceptor - donor);
                std::array<double, 3> const dh = (std::array<double, 3>) (*h - donor);
                double angle_adh = geom::angle<3>(da, dh);

                std::array<double, 3> const ha = (std::array<double, 3>) (acceptor - *h);
                std::array<double, 3> const hd = (std::array<double, 3>) (donor - *h);
                double angle_ahd = geom::angle<3>(ha, hd);

                if (angle_adh <= cfg::params::hbond_angle) // 63
                {
                    net.find(acceptor.res(), donor.res()).push(*new bond::hydrogen(params, acceptor, donor, h, angle_ahd));
                    return true;
                }
            }
        }
    }

    return false;
}

bool bond::vdw::test(network& net, rin::parameters const& params, chemical_entity::atom const& a, chemical_entity::atom const& b)
{
    if (a.res().satisfies_minimum_separation(b.res()) &&
        a.distance(b) - (a.vdw_radius() + b.vdw_radius()) <= params.surface_dist_vdw())
    {
        // FIXME we take the bonds two times
        auto& pb = net.find(a.res(), b.res());
        if (!pb.has_vdw())
            pb.push(*new bond::vdw(params, a, b));
        return true;
    }

    return false;
}


bool bond::ionic::test(network& net, rin::parameters const& params, chemical_entity::ionic_group const& negative, chemical_entity::ionic_group const& positive)
{
    if (negative.res().satisfies_minimum_separation(positive.res()) && negative.charge() == -positive.charge())
    {
        net.find(negative.res(), positive.res()).push(*new bond::ionic(params, negative, positive));
        return true;
    }

    return false;
}


bool bond::pication::test(network& net, rin::parameters const& params, chemical_entity::atom const& cation, chemical_entity::ring const& ring)
{
    if (ring.res().satisfies_minimum_separation(cation.res(), params.sequence_separation()))
    {
        double theta = 90 - geom::d_angle<3>(ring.normal(), (std::array<double, 3>) (ring - cation));

        if (theta >= cfg::params::pication_angle) // 45
        {
            net.find(ring.res(), cation.res()).push(*new bond::pication(params, ring, cation, theta));
            return true;
        }
    }

    return false;
}

bool bond::pipistack::test(network& net, rin::parameters const& params, chemical_entity::ring const& a, chemical_entity::ring const& b)
{
    double nc1 = a.angle_between_normal_and_centres_joining(b);
    double nc2 = b.angle_between_normal_and_centres_joining(a);
    double nn = a.angle_between_normals(b);
    double mn = a.closest_distance_between_atoms(b);

    if (a.res().satisfies_minimum_separation(b.res()) &&
        (0 <= nn && nn <= cfg::params::pipistack_normal_normal_angle_range) &&
        ((0 <= nc1 && nc1 <= cfg::params::pipistack_normal_centre_angle_range) ||
         (0 <= nc2 && nc2 <= cfg::params::pipistack_normal_centre_angle_range)) &&
        mn <= cfg::params::max_pipi_atom_atom_distance)
    {
        auto& pb = net.find(a.res(), b.res());
        if (!pb.has_pipi())
            pb.push(*new pipistack(params, a, b, nn));

        return true;
    }

    return false;
}

bool bond::generico::test(network& net, rin::parameters const& params, chemical_entity::atom const& a, chemical_entity::atom const& b)
{
    if (a.res().satisfies_minimum_separation(b.res()))
    {
        auto& pb = net.find(a.res(), b.res());
        if (!pb.has_backbone())
            pb.push(*new bond::generico(params, a.res(), b.res()));

        return true;
    }

    return false;
}