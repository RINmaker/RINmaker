#include "ns_bond.h"

#include <utility>
#include <map>
#include <set>

using chemical_entity::aminoacid, chemical_entity::atom, chemical_entity::ring, chemical_entity::ionic_group;
using rin::parameters;

using std::string, std::pair, std::make_pair, std::array, std::set;

// vv fixing "one/multiple" bug vv

template<typename Bond>
struct get;

template<>
struct get<bond::ss>
{
    static constexpr auto source_id = [](bond::ss const& bond)
        { return bond.get_source_id(); };

    static constexpr auto target_id = [](bond::ss const& bond)
    { return bond.get_target_id(); };
};

template<>
struct get<bond::ionic>
{
    static constexpr auto source_id = [](bond::ionic const& bond)
        { return bond.get_source_positive().get_residue().get_id(); };

    static constexpr auto target_id = [](bond::ionic const& bond)
    { return bond.get_target_negative().get_residue().get_id(); };
};

template<>
struct get<bond::hydrogen>
{
    static constexpr auto source_id = [](bond::hydrogen const& bond)
        { return bond.get_source_atom().get_residue().get_id(); };

    static constexpr auto target_id = [](bond::hydrogen const& bond)
    { return bond.get_target_atom().get_residue().get_id(); };
};

template<>
struct get<bond::pipistack>
{
    static constexpr auto source_id = [](bond::pipistack const& bond)
        { return bond.get_source_ring().get_residue().get_id(); };

    static constexpr auto target_id = [](bond::pipistack const& bond)
    { return bond.get_target_ring().get_residue().get_id(); };
};

template<>
struct get<bond::pication>
{
    static constexpr auto source_id = [](bond::pication const& bond)
        { return bond.get_source_ring().get_residue().get_id(); };

    static constexpr auto target_id = [](bond::pication const& bond)
    { return bond.get_target_cation().get_residue().get_id(); };
};

template<>
struct get<bond::vdw>
{
    static constexpr auto source_id = [](bond::vdw const& bond)
        { return bond.get_source_atom().get_residue().get_id(); };

    static constexpr auto target_id = [](bond::vdw const& bond)
    { return bond.get_target_atom().get_residue().get_id(); };
};

template<>
struct get<bond::generic_bond>
{
    static constexpr auto source_id = [](bond::generic_bond const& bond)
        { return bond.get_source().get_id(); };

    static constexpr auto target_id = [](bond::generic_bond const& bond)
        { return bond.get_target().get_id(); };
};

static constexpr auto lexicord = [](std::string const& a, std::string const& b)
{ return a < b ? a + b : b + a; };

template<typename Bond>
static std::string get_pair_id(Bond const& bond)
{
    using get = ::get<Bond>;
    auto s = get::source_id(bond);
    auto t = get::target_id(bond);
    return lexicord(s, t);
}

std::string bond::ss::get_pair_id() const
{
    return ::get_pair_id(*this);
}

std::string bond::vdw::get_pair_id() const
{
    return ::get_pair_id(*this);
}

std::string bond::hydrogen::get_pair_id() const
{
    return ::get_pair_id(*this);
}

std::string bond::generic_bond::get_pair_id() const
{
    return ::get_pair_id(*this);
}

std::string bond::pipistack::get_pair_id() const
{
    return ::get_pair_id(*this);
}

std::string bond::pication::get_pair_id() const
{
    return ::get_pair_id(*this);
}

std::string bond::ionic::get_pair_id() const
{
    return ::get_pair_id(*this);
}

// ^^ fixing "one/multiple" bug ^^

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

template<typename Entity>
bool operator<(Entity const& a, Entity const& b)
{ return a.get_residue().get_id() < b.get_residue().get_id(); }

generic_bond::generic_bond(atom const& a, atom const& b, double distance) :
    base(distance, 0), // todo: energy should be != 0 (?)
    _source(a < b ? a : b),
    _target(a < b ? b : a)
{}

chemical_entity::aminoacid generic_bond::get_source() const
{ return _source.get_residue(); }

chemical_entity::aminoacid generic_bond::get_target() const
{ return _target.get_residue(); }

string generic_bond::get_interaction() const
{ return "GENERIC:" + _source.get_name(); }

//Returns a pair of Sigmaij Epsilonij
pair<double, double> hydrogen::getSigmaEpsilon(atom const& donor, atom const& acceptor)
{
    std::function<bool(string const&, int, string const&, int)> compare =
        [&](string const& donor_element, int donor_charge, string const& acceptor_element, int acceptor_charge)
        {
            return donor.get_symbol() == donor_element
                && donor.get_charge() == donor_charge
                && acceptor.get_symbol() == acceptor_element
                && acceptor.get_charge() == acceptor_charge;
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

atom const& hydrogen::get_acceptor() const
{ return _acceptor; }

atom const& hydrogen::get_donor() const
{ return _donor; }

atom const& hydrogen::get_hydrogen_atom() const
{ return _hydrogen; }

double hydrogen::get_angle() const
{ return _angle; }

string hydrogen::get_interaction() const
{
    string donorChain = _donor.is_main_chain() ? "MC" : "SC";
    string acceptorChain = _acceptor.is_main_chain() ? "MC" : "SC";

    return "HBOND:" + acceptorChain + "_" + donorChain;
}

ionic::ionic(ionic_group const& negative, ionic_group const& positive) :
    base(negative.distance(positive), (constant::ion_ion_k * positive.get_ionion_energy_q() * negative.get_ionion_energy_q() / (negative.distance(positive)))),
    _negative(negative),
    _positive(positive)
{}

string ionic::get_interaction() const
{ return "IONIC:SC_SC"; }

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

double pication::get_angle() const
{ return _angle; }

string pication::get_interaction() const
{ return "PICATION:SC_SC"; }

double pipistack::energy(double angle)
{
    double cos_part = cos(1. / (angle + 10.));
    return constant::pipi_a + (constant::pipi_b * angle) + (constant::pipi_c * angle * cos_part);
}

pipistack::pipistack(ring const& a, ring const& b, double angle) :
    base(a.distance(b),energy(angle)),
    _source_ring(a < b ? a : b),
    _target_ring(a < b ? b : a),
    _angle(angle)
{}

double pipistack::get_angle() const
{ return _angle; }

string pipistack::get_interaction() const
{ return "PIPISTACK:SC_SC"; }

ss::ss(gemmi::Connection const& connection) :
    base(connection.reported_distance, 167),
    _source_seq{connection.partner1.res_id.seqid.num.value},
    _source_name{connection.partner1.res_id.name},
    _source_chain{connection.partner1.chain_name},

    _target_seq{connection.partner2.res_id.seqid.num.value},
    _target_name{connection.partner2.res_id.name},
    _target_chain{connection.partner2.chain_name}
{}

string ss::get_source_id() const
{ return _source_chain + ":" + std::to_string(_source_seq) + ":_:" + _source_name; }

string ss::get_target_id() const
{ return _target_chain + ":" + std::to_string(_target_seq) + ":_:" + _target_name; }

string ss::get_interaction() const
{ return "SSBOND:SC_SC"; } // TODO config

string ss::get_id() const
{ return get_id_simple(); }

double vdw::energy(atom const& source_atom, atom const& target_atom)
{
    double* source_opts = get_vdw_opsl_values(source_atom.get_residue().get_name(), source_atom.get_name(), source_atom.get_symbol());
    double* target_opts = get_vdw_opsl_values(target_atom.get_residue().get_name(), target_atom.get_name(), target_atom.get_symbol());

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
    base(a.distance(b), energy(a < b ? a : b, a < b ? b : a)),
    _source_atom(a < b ? a : b),
    _target_atom(a < b ? b : a)
{}

string vdw::get_interaction() const
{
    string sourceChain = _source_atom.get_name() == "C" || _source_atom.get_name() == "S" ? "MC" : "SC";
    string targetChain = _target_atom.get_name() == "C" || _target_atom.get_name() == "S" ? "MC" : "SC";

    return "VDW:" + sourceChain + "_" + targetChain;
}

std::shared_ptr<hydrogen const> hydrogen::test(parameters const& params, atom const& acceptor, atom const& donor)
{
    if (acceptor.get_residue().satisfies_minimum_sequence_separation(donor.get_residue()))
    {
        if (!(acceptor.get_residue() == donor.get_residue()))
        {
            auto hydrogens = donor.get_attached_hydrogens();
            for (auto const& h: hydrogens)
            {
                auto const da = (array<double, 3>) (acceptor - donor);
                auto const dh = (array<double, 3>) (h - donor);
                double angle_adh = geom::angle<3>(da, dh);

                auto const ha = (array<double, 3>) (acceptor - h);
                auto const hd = (array<double, 3>) (donor - h);
                double angle_ahd = geom::angle<3>(ha, hd);

                if (angle_adh <= params.hbond_angle()) // 63
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
        get_source_atom().get_name() +
        ":" +
        get_target_atom().get_name();
}

std::shared_ptr<vdw const> vdw::test(parameters const& params, atom const& a, atom const& b)
{
    if (a.get_residue().satisfies_minimum_sequence_separation(b.get_residue()) && a.distance(b) - (a.get_vdw_radius() + b.get_vdw_radius()) <= params.surface_dist_vdw())
        return std::make_shared<vdw const>(a, b);
    return nullptr;
}

string vdw::get_id() const
{
    /*
    return get_id_simple() +
        ":" +
        get_source_atom().get_name() +
        ":" +
        get_target_atom().get_name();
    */
    return get_id_simple() +
        ":" +
        get_source_atom().get_name() +
        ":" +
        std::to_string(get_source_atom().get_atom_number()) +
        ":" +
        get_target_atom().get_name() +
        ":" +
        std::to_string(get_target_atom().get_atom_number());
}


std::shared_ptr<ionic const> ionic::test(parameters const& params, ionic_group const& negative, ionic_group const& positive)
{
    if (negative.get_residue().satisfies_minimum_sequence_separation(positive.get_residue()) && negative.get_charge() == -positive.get_charge())
        return std::make_shared<ionic const>(negative, positive);

    return nullptr;
}

string ionic::get_id() const
{
    return get_id_simple() +
        ":" +
        get_source_positive().get_name() +
        ":" +
        get_target_negative().get_name();
}

std::shared_ptr<pication const> pication::test(parameters const& params, atom const& cation, ring const& ring)
{
    if (ring.get_residue().satisfies_minimum_sequence_separation(cation.get_residue(), params.sequence_separation()))
    {
        double theta = 90 - geom::d_angle<3>(ring.get_normal(), (array<double, 3>) (ring - cation));
        if (theta >= params.pication_angle()) // 45
            return std::make_shared<pication const>(ring, cation, theta);
    }

    return nullptr;
}

string pication::get_id() const
{
    return get_id_simple() +
        ":" +
        get_source_ring().get_name() +
        ":" +
        get_target_cation().get_name();
}

std::shared_ptr<pipistack const> pipistack::test(parameters const& params, ring const& a, ring const& b)
{
    double nc1 = a.get_angle_between_normal_and_centers_joining(b);
    double nc2 = b.get_angle_between_normal_and_centers_joining(a);
    double nn = a.get_angle_between_normals(b);
    double mn = a.get_distance_between_closest_atoms(b);

    if (a.get_residue().satisfies_minimum_sequence_separation(b.get_residue()) &&
        (0 <= nn && nn <= params.pipistack_normal_normal_angle_range()) &&
        ((0 <= nc1 && nc1 <= params.pipistack_normal_centre_angle_range()) ||
         (0 <= nc2 && nc2 <= params.pipistack_normal_centre_angle_range())) &&
        mn <= cfg::params::max_pipi_atom_atom_distance)
    { return std::make_shared<pipistack const>(a, b, nn); }

    return nullptr;
}

string pipistack::get_id() const
{
    return get_id_simple() +
        ":" +
        get_source_ring().get_name() +
        ":" +
        get_target_ring().get_name();
}

string generic_bond::get_id() const
{ return get_id_simple(); }

string generic_bond::get_id_simple() const
{
    return "GENERIC:" +
        get_source().get_id() +
        ":" +
        get_target().get_id();
}

string pipistack::get_id_simple() const
{
    return "PIPISTACK:" +
        get_source_ring().get_residue().get_id() +
        ":" +
        get_target_ring().get_residue().get_id();
}

string pication::get_id_simple() const
{
    return "PICATION:" +
        get_source_ring().get_residue().get_id() +
        ":" +
        get_target_cation().get_residue().get_id();
}

string hydrogen::get_id_simple() const
{
    return "HYDROGEN:" +
        get_source_atom().get_residue().get_id() +
        ":" +
        get_target_atom().get_residue().get_id() +
        ":" +
        get_hydrogen_atom().get_name();
}

string vdw::get_id_simple() const
{
    return "VDW:" +
        get_source_atom().get_residue().get_id() +
        ":" +
        get_target_atom().get_residue().get_id();
}

string ionic::get_id_simple() const
{
    return "IONIC:" +
        get_source_positive().get_residue().get_id() +
        ":" +
        get_target_negative().get_residue().get_id();
}

string ss::get_id_simple() const
{ return "SS:" + get_source_id() + ":" + get_target_id(); }

std::shared_ptr<hydrophobic const> hydrophobic::test(rin::parameters const&, chemical_entity::atom const& a, chemical_entity::atom const& b)
{
    static set<string> const names = {"ILE", "LEU", "VAL", "MET", "PHE", "ALA", "TRP", "CYS", "GLY"};
    auto res_a = a.get_residue();
    auto res_b = b.get_residue();

    if (res_a.satisfies_minimum_sequence_separation(res_b)
        && names.find(res_a.get_name()) != names.end()
        && names.find(res_b.get_name()) != names.end())
    { return std::make_shared<hydrophobic>(a, b); }

    return nullptr;
}

hydrophobic::hydrophobic(chemical_entity::atom const& c1, chemical_entity::atom const& c2) :
    generic_bond(c1, c2, geom::distance((std::array<double, 3>) c1, (std::array<double, 3>) c2))
{
    // Tab. 4 p. 13
    static std::map<string, double> const alpha = {
        {"ILE", 95.2},
        {"LEU", 94.5},
        {"VAL", 81.5},
        {"MET", 102.1},
        {"PHE", 122.9},
        {"ALA", 55.9},
        {"TRP", 157.8},
        {"CYS", 77}, // computed
        {"GLY", 44.3},
    };

    static constexpr auto get_alpha = [](auto const& atom)
        {
            auto it = alpha.find(atom.get_residue().get_name());
            if (it == alpha.end())
            {
                // At this stage it is necessary to "know" about the error (in theory) but it will never be thrown.
                // If it happens, fix immediately.
                // TODO: even if it does not happen, think of a better way.
                throw std::runtime_error("residue name out of mapping's domain");
            }
            return it->second;
        };

    auto const a1 = get_alpha(c1);
    auto const a2 = get_alpha(c2);

    static constexpr auto h = 6.6260700e-34;
    static constexpr auto n = 1.;
    static constexpr auto ni = 4e13 * n;
    static constexpr auto e0 = 8.854e-12;
    static constexpr auto pi = 3.1415927;

    _energy = -3*h*ni*a1*a2/(4*pow(4*pi*e0,2)*pow(_length,6));
}

std::string hydrophobic::get_interaction() const
{
    return "HYDROPHOBIC";
}

std::shared_ptr<contact const> contact::test(parameters const& params, atom const& a, atom const& b)
{
    if (a.get_residue().satisfies_minimum_sequence_separation(b.get_residue()))
        return std::make_shared<contact const>(a, b);
    return nullptr;
}

contact::contact(const atom& a, const atom& b) :
    generic_bond(a, b, geom::distance(a.get_residue().get_position(), b.get_residue().get_position()))
{}
