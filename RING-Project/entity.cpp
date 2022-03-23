#include "entity.h"
#include "prelude.h"
#include "geom.h"
#include "config.h"
#include "records.h"
#include "energy.h"

entities::aminoacid::aminoacid(std::vector<records::atom> const& records) : kdpoint<3>({0, 0, 0})
{
    /* TODO throw exception. It cannot happen.
    if (records.empty())
    {
    }
    */

    records::atom const& first = records.front();
    _name = first.res_name();
    _sequence_number = first.res_seq();
    _chain_id = first.chain_id();

    _id = _chain_id + ":" + std::to_string(_sequence_number) + ":_:" + _name;

    // discover if this has 0, 1 or 2 aromatic rings
    int n_of_rings = 0;
    std::vector<atom const*> ring_1, ring_2;
    std::vector<std::string> patterns_1, patterns_2;

    if (_name == "HIS")
    {
        // HIS has a 5-atoms ring
        patterns_1 = {"CD", "CE", "CG", "ND", "NE"};
        n_of_rings = 1;
    }
    else if (_name == "PHE" || _name == "TYR")
    {
        // PHE, TYR have a 6-atoms ring
        patterns_1 = {"CD1", "CD2", "CE1", "CE2", "CG", "CZ"};
        n_of_rings = 1;
    }
    else if (_name == "TRP")
    {
        // TRP has both a 6-atoms ring and a 5-atoms ring
        patterns_1 = {"CD2", "CE2", "CE3", "CH", "CZ2", "CZ3"};
        patterns_2 = {"CD2", "CE2", "CD1", "CG", "NE"};

        n_of_rings = 2;
    }

    // TODO are they mutually exclusive? should be addressed
    std::vector<atom const*> positive, negative;

    std::array<double, 3> centroid({0, 0, 0});
    double mass = 0;
    for (auto const& record: records)
    {
        atom const* a = new entities::atom(record, *this);
        _atoms.push_back(a);

        // sum up centroid weights
        double const m = a->mass();
        for (size_t i = 0; i < 3; ++i)
            centroid[i] += (*a)[i] * m;

        // sum up the total mass
        mass += m;

        if (a->name() == "CA")
        {
            // try to get alpha carbon (there should be only one)
            _alpha_carbon = a;
        }
        else if (a->name() == "CB")
        {
            // try to get alpha carbon (there should be only one, or zero)
            _beta_carbon = a;
        }

        // try to build rings

        if (n_of_rings >= 1 && prelude::match(a->name(), patterns_1))
        {
            ring_1.push_back(a);
        }
        if (n_of_rings == 2 && prelude::match(a->name(), patterns_2))
        {
            ring_2.push_back(a);
        }

        // try to build ionic groups
        if (a->is_in_a_positive_ionic_group())
        {
            positive.push_back(a);
        }
        else if (a->is_in_a_negative_ionic_group())
        {
            negative.push_back(a);
        }
    }

    for (size_t i = 0; i < 3; ++i)
        centroid[i] /= mass;

    _position = centroid;

    if (!ring_1.empty())
    {
        _primary_ring = new ring(ring_1, *this);
    }

    if (!ring_2.empty())
    {
        _secondary_ring = new ring(ring_2, *this);
    }
    if (!positive.empty())
    {
        _positive_ionic_group = new ionic_group(positive, 1, *this);
    }
    if (!negative.empty())
    {
        _negative_ionic_group = new ionic_group(negative, -1, *this);
    }
}

entities::aminoacid::~aminoacid()
{
    for (auto* a: _atoms)
        delete a;

    delete _primary_ring;
    delete _secondary_ring;

    delete _positive_ionic_group;
    delete _negative_ionic_group;

    delete _secondary_structure;
}

double entities::atom::mass() const
{
    std::string element = symbol();
    if (element == "H")
    {
        return 1.008;
    }
    else if (element == "C")
    {
        return 12.011;
    }
    else if (element == "N")
    {
        return 14.007;
    }
    else if (element == "O")
    {
        return 15.994;
    }
    else if (element == "S")
    {
        return 32.065;
    }

    // this point should be never reached (it depends on the guarantees of the contents of PDB files)
    // should it happen, we count this as "hydrogen"
    return 1;
}

double entities::atom::vdw_radius() const
{
    std::string element = symbol();
    if (element == "S")
    {
        return 1.89;
    }
    else if (element == "C")
    {
        return 1.77;
    }
    else if (element == "O")
    {
        return 1.55;
    }
    else if (element == "N")
    {
        return 1.60;
    }

    // chimicamente non ci aspettiamo di arrivare qui
    // dovesse succedere lo contiamo come "non-vdw"
    return 0;
}

bool entities::atom::is_a_cation() const
{
    std::string res_name = res().name();

    return (res_name == "LYS" && name() == "NZ") ||
           (res_name == "ARG" && name() == "NH") ||
           (res_name == "HIS" && name() == "ND1");
}

bool entities::atom::is_in_a_positive_ionic_group() const
{
    std::string res_name = res().name();

    if (res_name == "HIS")
    {
        return name() == "CG" || name() == "CD2" || name() == "CE1" || name() == "ND1" || name() == "NE2";
    }
    else if (res_name == "ARG")
    {
        return name() == "CZ" || name() == "NH2";
    }
    else if (res_name == "LYS")
    {
        return name() == "NZ";
    }
    return false;
}

bool entities::atom::is_in_a_negative_ionic_group() const
{
    std::string res_name = res().name();
    std::string n = name();

    if (res_name == "GLU")
    {
        return n == "CD" || n == "OE1" || n == "OE2";
    }

    else if (res_name == "ASP")
    {
        return n == "CG" || n == "OD1" || n == "OD2";
    }

    return false;
}

bool entities::atom::is_a_hydrogen_donor() const
{
    std::string res_name = res().name();
    std::string n = name();
    return
            (res_name == "ARG" && (n == "NH1" || n == "NH2" || n == "NE")) ||
            (res_name == "ASN" && n == "ND2") ||
            (res_name == "GLN" && n == "NE2") ||
            (res_name == "HIS" && (n == "NE2" || n == "ND1")) ||
            (res_name == "LYS" && n == "NZ") ||
            (res_name == "SER" && n == "OG") ||
            (res_name == "THR" && n == "OG1") ||
            (res_name == "TRP" && n == "NE1") ||
            (res_name == "TYR" && n == "OH") ||
            (res_name == "CYS" && n == "SG") ||
            n == "NH" || n == "N";
}
int entities::atom::how_many_hydrogen_can_donate() const
{
    if (is_a_hydrogen_donor())
    {
        std::string res_name = res().name();
        std::string n = name();
        if ((res_name == "ARG" && (n == "NH1" || n == "NH2")) ||
            (res_name == "ASN" && n == "ND2") ||
            (res_name == "GLN" && n == "NE2"))
            return 2;
        else if ((res_name == "LYS" && n == "NZ"))
            return 3;
        else
            return 1;
    }
    else
    {
        return 0;
    }
}

bool entities::atom::is_a_hydrogen_acceptor() const
{
    std::string res_name = res().name();
    std::string n = name();

    return
            (res_name == "ASN" && n == "OD1") ||
            (res_name == "ASP" && (n == "OD1" || n == "OD2")) ||
            (res_name == "GLN" && n == "OE1") ||
            (res_name == "GLU" && (n == "OE1" || n == "OE2")) ||
            (res_name == "HIS" && (n == "ND1" || n == "NE2")) ||
            (res_name == "SER" && n == "OG") ||
            (res_name == "THR" && n == "OG1") ||
            (res_name == "TYR" && n == "OH") ||
            (res_name == "MET" && n == "SD") ||
            n == "C" || n == "O";
}
int entities::atom::how_many_hydrogen_can_accept() const
{
    if (is_a_hydrogen_acceptor())
    {
        std::string res_name = res().name();
        std::string n = name();
        if ((res_name == "ASN" && n == "OD1") ||
            (res_name == "ASP" && (n == "OD1" || n == "OD2")) ||
            (res_name == "GLN" && n == "OE1") ||
            (res_name == "GLU" && (n == "OE1" || n == "OE2")) ||
            (res_name == "SER" && n == "OG") ||
            (res_name == "THR" && n == "OG1"))
            return 2;
        else
            return 1;
    }
    else
    {
        return 0;
    }
}

bool entities::atom::is_a_vdw_candidate() const
{
    std::string res_name = res().name();
    std::string n = name();
    std::string en = symbol();

    double* opslVdwValue = get_vdw_opsl_values(res_name, n, en);
    return !(opslVdwValue[0] == 0. &&
        opslVdwValue[1] == 0. &&
        opslVdwValue[2] == 0.);

    /*return
            (res_name == "GLN" && (n == "NE1" || n == "OE1")) ||
            (res_name == "ASN" && (n == "ND2" || n == "OD1")) ||
            en == "C" || en == "S";*/
}

std::vector<entities::atom const*> entities::atom::attached_hydrogens() const
{
    std::vector<atom const*> hydrogens;
    std::string const hydrogen_name_pattern = "H" + name().substr(1, name().size() - 1);
    for (atom const* a: res().atoms())
    {
        if (a->is_a_hydrogen() && prelude::match(a->name(), hydrogen_name_pattern))
        {
            hydrogens.push_back(a);
        }
    }

    return hydrogens;
}

entities::ring::ring(std::vector<atom const*> const& atoms, aminoacid const& res)
        : kdpoint<3>({0, 0, 0}), component(res), _atoms(atoms)
{
    // TODO assert atoms.size >= 3
    std::array<double, 3> centroid({0, 0, 0});
    double mass = 0;
    double sum_radii = 0;
    for (auto* a: atoms)
    {
        double m = a->mass();
        mass += m;
        sum_radii += distance(*a);

        for (size_t i = 0; i < 3; ++i)
            centroid[i] += (*a)[i] * m;
    }

    for (size_t i = 0; i < 3; ++i)
        centroid[i] /= mass;

    _position = centroid;
    _mean_radius = sum_radii / (double) atoms.size();

    // kudos to Giulio Marcolin for this shortcut
    // it only misses a SVD best-fit method no more than 1-2°, on average

    std::array<double, 3> const v = (std::array<double, 3>) ((*atoms[0]) - (*atoms[1]));
    std::array<double, 3> const w = (std::array<double, 3>) ((*atoms[2]) - (*atoms[1]));
    _normal = geom::cross(v, w);
}

entities::atom const& entities::ring::atom_closest_to(atom const& atom) const
{
    auto* closest_atom = _atoms[0];
    double min_distance = closest_atom->distance(atom);

    for (auto* a: _atoms)
    {
        const double distance = a->distance(atom);
        if (distance < min_distance)
        {
            min_distance = distance;
            closest_atom = a;
        }
    }

    return *closest_atom;
}

string getNameFromAtoms(std::vector<const entities::atom*> atoms)
{
    vector<string> atoms_name;
    for (auto i : atoms)
        atoms_name.push_back(i->name());

    sort(atoms_name.begin(), atoms_name.end());
    
    string out;
    string delimiter = ":";
    for (string& name : atoms_name)
    {
        if (!out.empty())
            out += delimiter;
        out += name;
    }
    return out;
}

string entities::ring::name() const
{
    return getNameFromAtoms(_atoms);
}

entities::ionic_group::ionic_group(std::vector<atom const*> const& atoms, int const& charge, aminoacid const& res)
        : kdpoint<3>({0, 0, 0}), component(res), _atoms(atoms), _charge(charge)
{
    // centroid is the centre of mass
    std::array<double, 3> centroid({0, 0, 0});
    double mass = 0;
    for (auto* a: atoms)
    {
        double m = a->mass();
        mass += m;

        for (size_t i = 0; i < 3; ++i)
            centroid[i] += (*a)[i] * m;
    }

    for (size_t i = 0; i < 3; ++i)
        centroid[i] /= mass;

    _position = centroid;
}

double entities::ionic_group::ionion_energy_q() const
{
    //                                    q  // * number of protons
    if (_res->name() == "LYS") return 0.640; // * 81;
    if (_res->name() == "ASP") return 0.380; // * 95;
    if (_res->name() == "HIS") return 0.380; // * 83;
    if (_res->name() == "ARG") return 0.260; // * 77;
    if (_res->name() == "GLU") return 0.635; // * 69;

    // TODO exception
    return 0;
}

string entities::ionic_group::name() const
{
    return getNameFromAtoms(_atoms);
}