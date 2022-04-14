#include "bonds.h"

#include <utility>

#include "chemical_entity.h"
#include "rin_network.h"

using chemical_entity::aminoacid, chemical_entity::atom, chemical_entity::ring, chemical_entity::ionic_group;
using rin::parameters;

using std::string, std::pair, std::make_pair, std::array;

using namespace bond;

base::base(double length, double energy) : _length(length), _energy(energy)
{}

double base::get_length() const
{ return _length; }

double base::get_energy() const
{ return _energy; }

bool base::operator<(base const& rhs) const
{ return _energy < rhs._energy || (_energy == rhs._energy && _length < rhs._length); }

bool base::operator>(base const& rhs) const
{ return rhs < *this; }

// TODO memoize?
template <typename Entity>
pair<Entity const*, Entity const*> sort_by_res_name(Entity const& a, Entity const& b)
{ return a.res().name() < b.res().name() ? make_pair(&a, &b) : make_pair(&b, &a); }

generico::generico(parameters const& params, atom const& a, atom const& b) :
    base(a.res().distance(b.res()), 0), // TODO
    _source(sort_by_res_name(a, b).first->res()),
    _target(sort_by_res_name(a, b).second->res()),
    _params(params)
{}

string generico::get_interaction() const
{
    string get_interaction = "GENERIC:";
    switch (_params.interaction_type())
    {
    case parameters::interaction_type_t::ALPHA_BACKBONE:
        get_interaction += "CA";
        break;

    case parameters::interaction_type_t::BETA_BACKBONE:
        get_interaction += "CB";
        break;

    case parameters::interaction_type_t::NONCOVALENT_BONDS:
        get_interaction += "CLOSEST";
        break;
    }

    return get_interaction;
}

string generico::get_type() const
{ return "generic"; } // TODO config

//Returns a pair of Sigmaij Epsilonij
pair<double, double> hydrogen::getSigmaEpsilon(atom const& donor, atom const& acceptor)
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

double hydrogen::energy(atom const& donor, atom const& acceptor, atom const* hydrogen)
{
    pair<double, double> sigmaEpsilon = getSigmaEpsilon(donor, acceptor);
    double sigma = sigmaEpsilon.first;
    double epsilon = sigmaEpsilon.second;
    double distance = hydrogen->distance(acceptor);

    double sigma_distance_12 = pow(sigma / distance, 12);
    double sigma_distance_10 = pow(sigma / distance, 10);

    return 4 * epsilon * (sigma_distance_12 - sigma_distance_10);
}

hydrogen::hydrogen(atom const& acceptor, atom const& donor, atom const* hydrogen, double angle) :
        base(acceptor.distance(donor), energy(donor, acceptor, hydrogen)),
        _acceptor(acceptor),
        _donor(donor),
        _hydrogen(hydrogen),
        _angle(angle)
{}

atom const& hydrogen::acceptor() const
{ return _acceptor; }

atom const& hydrogen::donor() const
{ return _donor; }

atom const* hydrogen::acceptor_ptr() const
{ return &_acceptor; }

atom const* hydrogen::donor_ptr() const
{ return &_donor; }

atom const* hydrogen::hydrogen_ptr() const
{ return _hydrogen; }

double hydrogen::get_angle() const
{ return _angle; }

string hydrogen::get_interaction() const
{
    string donorChain = _donor.is_main_chain() ? "MC" : "SC";
    string acceptorChain = _acceptor.is_main_chain() ? "MC" : "SC";

    return "HBOND:" + acceptorChain + "_" + donorChain;
}

string hydrogen::get_type() const
{ return "hydrogen"; } // TODO config

ionic::ionic(ionic_group const& negative, ionic_group const& positive) :
        base(
                negative.distance(positive),
                (constant::ion_ion_k * positive.ionion_energy_q() * negative.ionion_energy_q() / (negative.distance(positive)))),
        _negative(negative),
        _positive(positive)
{}

string ionic::get_interaction() const
{ return "IONIC:SC_SC"; }

string ionic::get_type() const
{ return "ionic"; }

pication::pication(ring const& ring, atom const& cation, double angle) :
        base(ring.distance(cation), 9.6),// TODO add formula
        _cation(cation),
        _ring(ring),
        _angle(angle)
{}

double pication::angle() const
{ return _angle; }

string pication::get_interaction() const
{ return "PICATION:SC_SC"; }

string pication::get_type() const
{ return "pication"; }

pipistack::pipistack(ring const& a, ring const& b, double angle) :
        base(a.distance(b), 9.6), // TODO add formula

        _source_ring(*sort_by_res_name(a, b).first),
        _target_ring(*sort_by_res_name(a, b).second),

        _angle(angle)
{}

double pipistack::angle() const
{ return _angle; }

string pipistack::get_interaction() const
{ return "PIPISTACK:SC_SC"; }

string pipistack::get_type() const
{ return "pipistack"; }

ss::ss(records::ss const& record)
        : base(record.length(), 167), // TODO config
          _source_name(record.name_1()),
          _target_name(record.name_2()),
          _source_chain(record.chain_id_1()),
          _target_chain(record.chain_id_2()),
          _source_seq(record.seq_num_1()),
          _target_seq(record.seq_num_2())
{}

string ss::source_id() const
{ return _source_chain + ":" + std::to_string(_source_seq) + ":_:" + _source_name; }

string ss::target_id() const
{ return _target_chain + ":" + std::to_string(_target_seq) + ":_:" + _target_name; }

string ss::get_interaction() const
{ return "SSBOND:SC_SC"; } // TODO config

string ss::get_id() const
{ return prelude::concat_lexicographically(source_id(), target_id()); }

string ss::get_type() const
{ return "ss"; }


double vdw::energy(atom const& source_atom, atom const& target_atom)
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

vdw::vdw(atom const& a, atom const& b) :
        base(
                a.distance(b),
                energy(
                        *sort_by_res_name(a, b).first,
                        *sort_by_res_name(a, b).second)),

        _source_atom(*sort_by_res_name(a, b).first),
        _target_atom(*sort_by_res_name(a, b).second)
{}

string vdw::get_interaction() const
{
    string sourceChain = _source_atom.name() == "C" || _source_atom.name() == "S" ? "MC" : "SC";
    string targetChain = _target_atom.name() == "C" || _target_atom.name() == "S" ? "MC" : "SC";

    return "VDW:" + sourceChain + "_" + targetChain;
}

string vdw::get_type() const
{ return "vdw"; }

bool hydrogen::test(network& net, parameters const& params, atom const& acceptor, atom const& donor)
{
    if (acceptor.res().satisfies_minimum_separation(donor.res()))
    {
        if (!(acceptor.res() == donor.res()))
        {
            auto hydrogens = donor.attached_hydrogens();
            for (auto* h: hydrogens)
            {
                array<double, 3> const da = (array<double, 3>) (acceptor - donor);
                array<double, 3> const dh = (array<double, 3>) (*h - donor);
                double angle_adh = geom::angle<3>(da, dh);

                array<double, 3> const ha = (array<double, 3>) (acceptor - *h);
                array<double, 3> const hd = (array<double, 3>) (donor - *h);
                double angle_ahd = geom::angle<3>(ha, hd);

                if (angle_adh <= cfg::params::hbond_angle) // 63
                {
                    net.find(acceptor.res(), donor.res()).push(*new bond::hydrogen(acceptor, donor, h, angle_ahd));
                    return true;
                }
            }
        }
    }

    return false;
}

string hydrogen::get_id() const
{
    return "HYDROGEN:" +
           source_atom().res().name() +
           ":" +
           source_atom().name() +
           ":" +
           target_atom().res().name() +
           ":" +
           target_atom().name();
}

bool vdw::test(network& net, parameters const& params, atom const& a, atom const& b)
{
    if (a.res().satisfies_minimum_separation(b.res()) &&
        a.distance(b) - (a.vdw_radius() + b.vdw_radius()) <= params.surface_dist_vdw())
    {
        // FIXME we take the bonds two times
        auto& pb = net.find(a.res(), b.res());
        if (!pb.has_vdw())
            pb.push(*new bond::vdw(a, b));
        return true;
    }

    return false;
}

string vdw::get_id() const
{
    return "VDW:" +
        source_atom().res().name() +
        ":" +
        source_atom().name() +
        ":" +
        target_atom().res().name() +
        ":" +
        target_atom().name();
}


bool ionic::test(network& net, parameters const& params, ionic_group const& negative, ionic_group const& positive)
{
    if (negative.res().satisfies_minimum_separation(positive.res()) && negative.charge() == -positive.charge())
    {
        net.find(negative.res(), positive.res()).push(*new bond::ionic(negative, positive));
        return true;
    }

    return false;
}

string ionic::get_id() const
{
    return "IONIC:" +
           source_positive().res().name() +
           ":" +
           source_positive().name() +
           ":" +
           target_negative().res().name() +
           ":" +
           target_negative().name();
}


bool pication::test(network& net, parameters const& params, atom const& cation, ring const& ring)
{
    if (ring.res().satisfies_minimum_separation(cation.res(), params.sequence_separation()))
    {
        double theta = 90 - geom::d_angle<3>(ring.normal(), (array<double, 3>) (ring - cation));

        if (theta >= cfg::params::pication_angle) // 45
        {
            net.find(ring.res(), cation.res()).push(*new bond::pication(ring, cation, theta));
            return true;
        }
    }

    return false;
}

string pication::get_id() const
{
    return "PICATION:" +
           source_ring().res().name() +
           ":" +
           source_ring().name() +
           ":" +
           target_cation().res().name() +
           ":" +
           target_cation().name();
}

bool pipistack::test(network& net, parameters const& params, ring const& a, ring const& b)
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
            pb.push(*new pipistack(a, b, nn));

        return true;
    }

    return false;
}

string pipistack::get_id() const
{
    return "PIPISTACK:" +
        source_ring().res().name() +
        ":" +
        source_ring().name() +
        ":" +
        target_ring().res().name() +
        ":" +
        target_ring().name();
}

bool generico::test(network& net, parameters const& params, atom const& a, atom const& b)
{
    if (a.res().satisfies_minimum_separation(b.res()))
    {
        auto& pb = net.find(a.res(), b.res());
        if (!pb.has_backbone())
            pb.push(*new bond::generico(params, a, b));

        return true;
    }

    return false;
}

string generico::get_id() const
{
    return "GENERICO:" +
           source().name() +
           ":" +
           target().name();
}
