#include "ns_bond.h"

#include <utility>

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
template<typename Entity>
pair<Entity const*, Entity const*> sort_by_res_id(Entity const& a, Entity const& b)
{ return a.get_residue().get_id() < b.get_residue().get_id() ? make_pair(&a, &b) : make_pair(&b, &a); }

generic_bond::generic_bond(parameters const& params, atom const& a, atom const& b) :
        base(geom::distance(a.get_residue().get_position(), b.get_residue().get_position()), 0), // TODO
        _source(*sort_by_res_id(a, b).first),
        _target(*sort_by_res_id(a, b).second)
{}

chemical_entity::aminoacid generic_bond::source() const
{ return _source.get_residue(); }

chemical_entity::aminoacid generic_bond::target() const
{ return _target.get_residue(); }

string generic_bond::get_interaction() const
{
    return "GENERIC:" + _source.name();
}

string generic_bond::get_type() const
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

    return std::make_pair(1.79, -4.25); //Default value (is valid for all MC_MC bond)

//    string donor_repr = donor.name() + "(" + std::to_string(donor.charge()) + ")";
//    string acceptor_repr = acceptor.name() + "(" + std::to_string(acceptor.charge()) + ")";
//    throw std::invalid_argument("hydrogen::getSigmaEpsilon: donor " + donor_repr + ", acceptor " + acceptor_repr + " unsupported");
}

double hydrogen::energy(atom const& donor, atom const& acceptor, atom const& hydrogen)
{
    pair<double, double> sigmaEpsilon = getSigmaEpsilon(donor, acceptor);
    double sigma = sigmaEpsilon.first;
    double epsilon = sigmaEpsilon.second;
    double distance = hydrogen.distance(acceptor);

    double sigma_distance_12 = pow(sigma / distance, 12);
    double sigma_distance_10 = pow(sigma / distance, 10);

    return 4 * epsilon * (sigma_distance_12 - sigma_distance_10);
}

hydrogen::hydrogen(atom const& acceptor, atom const& donor, atom const& hydrogen, double angle) :
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



atom const& hydrogen::hydrogen_atom() const
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
        base(negative.distance(positive), (constant::ion_ion_k * positive.ionion_energy_q() * negative.ionion_energy_q() / (negative.distance(positive)))),
        _negative(negative),
        _positive(positive)
{}

string ionic::get_interaction() const
{ return "IONIC:SC_SC"; }

string ionic::get_type() const
{ return "ionic"; }

double pication::getKappa(atom const& cation)
{
    string res_name = cation.get_residue().get_name();

    if (res_name == "LYS" || res_name == "HIS") return 1.00;
    if (res_name == "ARG") return 0.25;

    throw std::invalid_argument("pication::getKappa: cation res name " + res_name + " unsupported");
}

double pication::getAlpha(ring const& ring)
{
    string res_name = ring.get_residue().get_name();

    if (res_name == "PHE" || res_name == "TYR") return 190;
    if (res_name == "TRP") return 150;

    throw std::invalid_argument("pication::getAlpha: ring res name " + res_name + " unsupported");
}

double pication::energy(ring const& ring, atom const& cation)
{
    double distance = ring.distance(cation);
    double kappa = getKappa(cation);
    double alpha = getAlpha(ring);


    double energy = -(kappa*alpha)/pow(distance, 4);
    return energy;
}

pication::pication(ring const& ring, atom const& cation, double angle) :
        base(ring.distance(cation), energy(ring, cation)),
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

double pipistack::energy(double angle)
{
    double cos_part = cos(1. / (angle + 10.));
    return constant::pipi_a + (constant::pipi_b * angle) + (constant::pipi_c*angle * cos_part);
}

pipistack::pipistack(ring const& a, ring const& b, double angle) :
        base(a.distance(b),energy(angle)),

        _source_ring(*sort_by_res_id(a, b).first),
        _target_ring(*sort_by_res_id(a, b).second),

        _angle(angle)
{}

double pipistack::angle() const
{ return _angle; }

string pipistack::get_interaction() const
{ return "PIPISTACK:SC_SC"; }

string pipistack::get_type() const
{ return "pipistack"; }

string ss::source_id() const
{ return _source_chain + ":" + std::to_string(_source_seq) + ":_:" + _source_name; }

string ss::target_id() const
{ return _target_chain + ":" + std::to_string(_target_seq) + ":_:" + _target_name; }

string ss::get_interaction() const
{ return "SSBOND:SC_SC"; } // TODO config

string ss::get_id() const
{ return get_id_simple(); }

string ss::get_type() const
{ return "ss"; }


double vdw::energy(atom const& source_atom, atom const& target_atom)
{
    double* source_opts = get_vdw_opsl_values(source_atom.get_residue().get_name(), source_atom.name(), source_atom.symbol());
    double* target_opts = get_vdw_opsl_values(target_atom.get_residue().get_name(), target_atom.name(), target_atom.symbol());

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
        base(a.distance(b), energy(*sort_by_res_id(a, b).first, *sort_by_res_id(a, b).second)),

        _source_atom(*sort_by_res_id(a, b).first),
        _target_atom(*sort_by_res_id(a, b).second)
{}

string vdw::get_interaction() const
{
    string sourceChain = _source_atom.name() == "C" || _source_atom.name() == "S" ? "MC" : "SC";
    string targetChain = _target_atom.name() == "C" || _target_atom.name() == "S" ? "MC" : "SC";

    return "VDW:" + sourceChain + "_" + targetChain;
}

string vdw::get_type() const
{ return "vdw"; }

std::shared_ptr<hydrogen const> hydrogen::test(parameters const& params, atom const& acceptor, atom const& donor)
{
    if (acceptor.get_residue().satisfies_minimum_sequence_separation(donor.get_residue()))
    {
        if (!(acceptor.get_residue() == donor.get_residue()))
        {
            auto hydrogens = donor.attached_hydrogens();
            for (auto const& h: hydrogens)
            {
                array<double, 3> const da = (array<double, 3>) (acceptor - donor);
                array<double, 3> const dh = (array<double, 3>) (h - donor);
                double angle_adh = geom::angle<3>(da, dh);

                array<double, 3> const ha = (array<double, 3>) (acceptor - h);
                array<double, 3> const hd = (array<double, 3>) (donor - h);
                double angle_ahd = geom::angle<3>(ha, hd);

                if (angle_adh <= cfg::params::hbond_angle) // 63
                    return std::make_shared<hydrogen const>(acceptor, donor, h, angle_ahd);
            }
        }
    }

    return nullptr;
}

string hydrogen::get_id() const
{
    return get_id_simple() +
           ":" +
           source_atom().name() +
           ":" +
           target_atom().name();
}

std::shared_ptr<vdw const> vdw::test(parameters const& params, atom const& a, atom const& b)
{
    if (a.get_residue().satisfies_minimum_sequence_separation(b.get_residue()) && a.distance(b) - (a.vdw_radius() + b.vdw_radius()) <= params.surface_dist_vdw())
        return std::make_shared<vdw const>(a, b);
    return nullptr;
}

string vdw::get_id() const
{
    return get_id_simple() +
           ":" +
           source_atom().name() +
           ":" +
           target_atom().name();
}


std::shared_ptr<ionic const> ionic::test(parameters const& params, ionic_group const& negative, ionic_group const& positive)
{
    if (negative.get_residue().satisfies_minimum_sequence_separation(positive.get_residue()) && negative.charge() == -positive.charge())
        return std::make_shared<ionic const>(negative, positive);

    return nullptr;
}

string ionic::get_id() const
{
    return get_id_simple() +
           ":" +
           source_positive().name() +
           ":" +
           target_negative().name();
}

std::shared_ptr<pication const> pication::test(parameters const& params, atom const& cation, ring const& ring)
{
    if (ring.get_residue().satisfies_minimum_sequence_separation(cation.get_residue(), params.sequence_separation()))
    {
        double theta = 90 - geom::d_angle<3>(ring.get_normal(), (array<double, 3>) (ring - cation));
        if (theta >= cfg::params::pication_angle) // 45
            return std::make_shared<pication const>(ring, cation, theta);
    }

    return nullptr;
}

string pication::get_id() const
{
    return get_id_simple() +
           ":" +
        source_ring().get_name() +
           ":" +
           target_cation().name();
}

std::shared_ptr<pipistack const> pipistack::test(parameters const& params, ring const& a, ring const& b)
{
    double nc1 = a.get_angle_between_normal_and_centers_joining(b);
    double nc2 = b.get_angle_between_normal_and_centers_joining(a);
    double nn = a.get_angle_between_normals(b);
    double mn = a.get_distance_between_closest_atoms(b);

    if (a.get_residue().satisfies_minimum_sequence_separation(b.get_residue()) &&
        (0 <= nn && nn <= cfg::params::pipistack_normal_normal_angle_range) &&
        ((0 <= nc1 && nc1 <= cfg::params::pipistack_normal_centre_angle_range) ||
         (0 <= nc2 && nc2 <= cfg::params::pipistack_normal_centre_angle_range)) &&
        mn <= cfg::params::max_pipi_atom_atom_distance)
    { return std::make_shared<pipistack const>(a, b, nn); }

    return nullptr;
}

string pipistack::get_id() const
{
    return get_id_simple() +
           ":" +
        source_ring().get_name() +
           ":" +
        target_ring().get_name();
}

std::shared_ptr<generic_bond const> generic_bond::test(parameters const& params, atom const& a, atom const& b)
{
    if (a.get_residue().satisfies_minimum_sequence_separation(b.get_residue()))
        return std::make_shared<generic_bond const>(params, a, b);;
    return nullptr;
}

string generic_bond::get_id() const
{
    return get_id_simple();
}

string generic_bond::get_id_simple() const
{
    return "GENERICO:" +
        source().get_id() +
           ":" +
        target().get_id();
}

string pipistack::get_id_simple() const
{
    return "PIPISTACK:" +
        source_ring().get_residue().get_id() +
           ":" +
        target_ring().get_residue().get_id();
}

string pication::get_id_simple() const
{
    return "PICATION:" +
        source_ring().get_residue().get_id() +
           ":" +
        target_cation().get_residue().get_id();
}

string hydrogen::get_id_simple() const
{
    return "HYDROGEN:" +
        source_atom().get_residue().get_id() +
           ":" +
        target_atom().get_residue().get_id() +
           ":" +
            hydrogen_atom().name();
}

string vdw::get_id_simple() const
{
    return "VDW:" +
        source_atom().get_residue().get_id() +
           ":" +
        target_atom().get_residue().get_id();
}

string ionic::get_id_simple() const
{
    return "IONIC:" +
        source_positive().get_residue().get_id() +
           ":" +
        target_negative().get_residue().get_id();
}

string ss::get_id_simple() const
{
    return "SS:" + source_id() + ":" + target_id();
}