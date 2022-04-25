#include "chemical_entity.h"

#include <memory>

#include "energy.h"

#include "prelude.h"
#include "utils.h"

using std::vector, std::array, std::string, std::unique_ptr, std::make_unique, std::to_string, std::invalid_argument;

using chemical_entity::aminoacid, chemical_entity::atom, chemical_entity::ring, chemical_entity::ionic_group,
        structure::helix, structure::loop, structure::sheet_piece;

string getNameFromAtoms(vector<atom const*> const& atoms, string const& delimiter = ":")
{
    vector<string> atoms_name;
    for (auto i: atoms)
        atoms_name.push_back(i->name());

    sort(atoms_name.begin(), atoms_name.end());

    return joinStrings(atoms_name, delimiter);
}

struct aminoacid::impl final
{
public:
    vector<unique_ptr<atom const>> _atoms;

    unique_ptr<ring const> _primary_ring, _secondary_ring;
    unique_ptr<ionic_group const> _positive_ionic_group, _negative_ionic_group;

    unique_ptr<structure::base> _secondary_structure{make_unique<structure::base>()};

    atom const* _alpha_carbon = nullptr;
    atom const* _beta_carbon = nullptr;

    string _chain_id;

    string _name;

    int _sequence_number = 0;

    string _id;

    string _pdb_name;

    ~impl() = default;
};

vector<atom const*> aminoacid::atoms() const
{
    vector<atom const*> obs;
    obs.reserve(pimpl->_atoms.size());

    for (auto const& a_uptr: pimpl->_atoms)
        obs.push_back(a_uptr.get());

    return obs;
}

string const& aminoacid::pdb_name() const
{ return pimpl->_pdb_name; }

atom const* aminoacid::ca() const
{ return pimpl->_alpha_carbon; }

atom const* aminoacid::cb() const
{ return pimpl->_beta_carbon; }

ring const* aminoacid::primary_ring() const
{ return pimpl->_primary_ring.get(); }

ring const* aminoacid::secondary_ring() const
{ return pimpl->_secondary_ring.get(); }

ionic_group const* aminoacid::positive_ionic_group() const
{ return pimpl->_positive_ionic_group.get(); }

ionic_group const* aminoacid::negative_ionic_group() const
{ return pimpl->_negative_ionic_group.get(); }

string const& aminoacid::name() const
{ return pimpl->_name; }

string const& aminoacid::chain_id() const
{ return pimpl->_chain_id; }

string const& aminoacid::id() const
{ return pimpl->_id; }

int aminoacid::sequence_number() const
{ return pimpl->_sequence_number; }

bool aminoacid::operator==(aminoacid const& rhs) const
{ return pimpl->_id == rhs.pimpl->_id; }

bool aminoacid::operator!=(aminoacid const& rhs) const
{ return !(*this == rhs); }

bool aminoacid::satisfies_minimum_separation(aminoacid const& aa, int minimum_separation) const
{
    if (*this == aa)
    {
        return false;
    }

    if (pimpl->_chain_id != aa.pimpl->_chain_id)
    {
        return true;
    }

    return abs(pimpl->_sequence_number - aa.pimpl->_sequence_number) >= minimum_separation;
}

aminoacid::operator rin::node() const
{ return rin::node(*this); }

void aminoacid::make_secondary_structure()
{
    pimpl->_secondary_structure = std::make_unique<structure::loop>();
}

void aminoacid::make_secondary_structure(records::helix const& record)
{
    pimpl->_secondary_structure = std::make_unique<structure::helix>(record, *this);
}

void aminoacid::make_secondary_structure(const records::sheet_piece& record)
{
    pimpl->_secondary_structure = std::make_unique<structure::sheet_piece>(record, *this);
}

string aminoacid::secondary_structure_id() const
{
    return pimpl->_secondary_structure->pretty();
}

array<double, 3> centre_of_mass(vector<atom const*> const& atoms)
{
    double mass{0.0};
    array<double, 3> centroid{0.0, 0.0, 0.0};
    for (auto a_ptr: atoms)
    {
        double const m = a_ptr->mass();
        for (size_t i = 0; i < 3; ++i)
            centroid[i] += (*a_ptr)[i] * m;
        mass += m;
    }
    for (size_t i = 0; i < 3; ++i)
        centroid[i] /= mass;

    return centroid;
}

aminoacid::aminoacid(vector<records::atom> const& records, string const& pdb_name) :
        kdpoint<3>({0, 0, 0}), pimpl{new impl()}
{
    auto assert_ring_correctness =
            [](string const& name, uint32_t line_number, vector<string> const& expected_atoms,
               vector<atom const*> const& found_atoms)
            {
                if (expected_atoms.size() != found_atoms.size())
                {
                    string expected_atoms_str = joinStrings(expected_atoms, ",");
                    string found_atoms_str = getNameFromAtoms(found_atoms, ",");
                    string exception_description =
                            "line number: " + std::to_string(line_number) + ", aminoacid: " + name +
                            " - expected aromatic ring: " + expected_atoms_str + ", found: " + found_atoms_str;
                    throw invalid_argument(exception_description);
                }
            };

    /* TODO throw exception. It cannot happen.
    if (records.empty())
    {
    }
    */

    auto const& first = records.front();
    pimpl->_name = first.res_name();
    pimpl->_sequence_number = first.res_seq();
    pimpl->_chain_id = first.chain_id();

    pimpl->_id = pimpl->_chain_id + ":" + to_string(pimpl->_sequence_number) + ":_:" + pimpl->_name;

    // discover if this has 0, 1 or 2 aromatic rings
    int n_of_rings = 0;
    vector<string> patterns_1, patterns_2;

    if (pimpl->_name == "HIS")
    {
        // HIS has a 5-atoms ring
        patterns_1 = {"CD2", "CE1", "CG", "ND1", "NE2"};
        n_of_rings = 1;
    }
    else if (pimpl->_name == "PHE" || pimpl->_name == "TYR")
    {
        // PHE, TYR have a 6-atoms ring
        patterns_1 = {"CD1", "CD2", "CE1", "CE2", "CG", "CZ"};
        n_of_rings = 1;
    }
    else if (pimpl->_name == "TRP")
    {
        // TRP has both a 6-atoms ring and a 5-atoms ring
        patterns_1 = {"CD2", "CE2", "CE3", "CH2", "CZ2", "CZ3"};
        patterns_2 = {"CD2", "CE2", "CD1", "CG", "NE1"};

        n_of_rings = 2;
    }

    // note: are they mutually exclusive? should be addressed
    vector<atom const*> positive, negative;
    vector<atom const*> ring_1, ring_2;

    for (auto const& record: records)
    {
        pimpl->_atoms.push_back(make_unique<atom const>(record, *this));
        auto const a = pimpl->_atoms.back().get();

        if (a->name() == "CA")
            pimpl->_alpha_carbon = a;

        else if (a->name() == "CB")
            pimpl->_beta_carbon = a;

        if (n_of_rings >= 1 && find(patterns_1.begin(), patterns_1.end(), a->name()) != patterns_1.end())
            ring_1.push_back(a);

        if (n_of_rings == 2 && find(patterns_2.begin(), patterns_2.end(), a->name()) != patterns_2.end())
            ring_2.push_back(a);

        if (a->is_in_a_positive_ionic_group())
            positive.push_back(a);

        else if (a->is_in_a_negative_ionic_group())
            negative.push_back(a);
    }

    _position = centre_of_mass(atoms());

    if (n_of_rings >= 1)
    {
        assert_ring_correctness(pimpl->_name, first.line_number(), patterns_1, ring_1);
        pimpl->_primary_ring = make_unique<ring const>(ring_1, *this);
    }
    if (n_of_rings == 2)
    {
        assert_ring_correctness(pimpl->_name, first.line_number(), patterns_2, ring_2);
        pimpl->_secondary_ring = make_unique<ring const>(ring_2, *this);
    }

    if (!positive.empty())
        pimpl->_positive_ionic_group = make_unique<ionic_group const>(positive, 1, *this);

    if (!negative.empty())
        pimpl->_negative_ionic_group = make_unique<ionic_group const>(negative, -1, *this);
}

aminoacid::~aminoacid()
{ delete pimpl; }

struct atom::impl final
{
public:
    records::atom _record;
};

atom::atom(records::atom const& record, aminoacid const& res) :
        kdpoint<3>({record.x(), record.y(), record.z()}), component(res), pimpl{new impl{record}}
{}

atom::~atom()
{ delete pimpl; }

string const& atom::name() const
{ return pimpl->_record.name(); }

string const& atom::symbol() const
{ return pimpl->_record.element_name(); }

double atom::temp_factor() const
{ return pimpl->_record.temp_factor(); }

int atom::charge() const
{ return pimpl->_record.charge(); }

bool atom::is_a_hydrogen() const
{ return symbol() == "H"; }

bool atom::is_main_chain() const
{
    return
            this->name() == "C" ||
            this->name() == "O" ||
            this->name() == "H" ||
            this->name() == "HA" ||
            this->name() == "N";
}

double atom::mass() const
{
    auto element = symbol();
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
    // TODO exception
    return 1;
}

double atom::vdw_radius() const
{
    string element = symbol();
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

    // TODO exception
    return 0;
}

bool atom::is_a_cation() const
{
    auto res_name = res().name();

    return (res_name == "LYS" && name() == "NZ") ||
           (res_name == "ARG" && name() == "NH2") ||
           (res_name == "HIS" && name() == "ND1");
}

bool atom::is_in_a_positive_ionic_group() const
{
    auto res_name = res().name();

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

bool chemical_entity::atom::is_in_a_negative_ionic_group() const
{
    auto res_name = res().name();
    auto n = name();

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

bool atom::is_a_hydrogen_donor() const
{
    auto res_name = res().name();
    auto n = name();
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

int atom::how_many_hydrogen_can_donate() const
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

bool atom::is_a_hydrogen_acceptor() const
{
    auto res_name = res().name();
    auto n = name();

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

int atom::how_many_hydrogen_can_accept() const
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

bool atom::is_a_vdw_candidate() const
{
    auto res_name = res().name();
    auto n = name();
    auto en = symbol();

    double* opslVdwValue = get_vdw_opsl_values(res_name, n, en);
    return !(opslVdwValue[0] == 0. &&
             opslVdwValue[1] == 0. &&
             opslVdwValue[2] == 0.);

    /*return
            (res_name == "GLN" && (n == "NE1" || n == "OE1")) ||
            (res_name == "ASN" && (n == "ND2" || n == "OD1")) ||
            en == "C" || en == "S";*/
}

vector<atom const*> atom::attached_hydrogens() const
{
    vector<atom const*> hydrogens;
    auto const hydrogen_name_pattern = "H" + name().substr(1, name().size() - 1);
    for (atom const* a: res().atoms())
    {
        if (a->is_a_hydrogen() && prelude::match(a->name(), hydrogen_name_pattern))
        {
            hydrogens.push_back(a);
        }
    }

    return hydrogens;
}

struct ring::impl final
{
public:
    vector<atom const*> _atoms;

    array<double, 3> _normal{};
    double _mean_radius;

};

ring::ring(vector<atom const*> const& atoms, aminoacid const& res) :
        kdpoint<3>({0, 0, 0}), component(res), pimpl{new impl()}
{
    if (atoms.size() < 3)
        throw invalid_argument("rings should have at least 3 atoms");

    pimpl->_atoms = atoms;

    double sum_radii = 0;
    for (auto* a: atoms)
        sum_radii += distance(*a);
    pimpl->_mean_radius = sum_radii / (double) atoms.size();

    _position = centre_of_mass(atoms);

    // kudos to Giulio Marcolin for the following shortcut
    // it only deviates from a SVD best-fit method no more than 1-2°, on average
    array<double, 3> const v = (array<double, 3>) ((*atoms[0]) - (*atoms[1]));
    array<double, 3> const w = (array<double, 3>) ((*atoms[2]) - (*atoms[1]));
    pimpl->_normal = geom::cross(v, w);
}

ring::~ring()
{ delete pimpl; }

array<double, 3> const& ring::normal() const
{ return pimpl->_normal; }

double ring::radius() const
{ return pimpl->_mean_radius; }

bool ring::is_a_pication_candidate() const
{
    string name = res().name();
    return name == "PHE" || name == "TYR" || (name == "TRP" && pimpl->_atoms.size() == 6);
}

double ring::closest_distance_between_atoms(ring const& other) const
{
    double minimum = pimpl->_atoms[0]->distance(*other.pimpl->_atoms[0]);
    for (auto* atom_1: pimpl->_atoms)
    {
        for (auto* atom_2: other.pimpl->_atoms)
        {
            double current = atom_1->distance(*atom_2);
            if (current < minimum)
            {
                minimum = current;
            }
        }
    }

    return minimum;
}

double ring::angle_between_normals(ring const& other) const
{ return geom::d_angle<3>(pimpl->_normal, other.pimpl->_normal); }

double ring::angle_between_normal_and_centres_joining(ring const& other) const
{
    std::array<double, 3> const centres_joining((std::array<double, 3>) (*this - other));
    return geom::d_angle<3>(pimpl->_normal, centres_joining);
}

atom const& ring::atom_closest_to(atom const& atom) const
{
    auto* closest_atom = pimpl->_atoms[0];
    double min_distance = closest_atom->distance(atom);

    for (auto* a: pimpl->_atoms)
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

string ring::name() const
{
    return getNameFromAtoms(pimpl->_atoms);
}

struct ionic_group::impl
{
public:
    vector<atom const*> const _atoms;
    int const _charge;
};

ionic_group::ionic_group(vector<atom const*> const& atoms, int const& charge, aminoacid const& res) :
        kdpoint<3>({0, 0, 0}), component(res), pimpl{new impl{atoms, charge}}
{ _position = centre_of_mass(atoms); }

ionic_group::~ionic_group()
{ delete pimpl; }

int ionic_group::charge() const
{ return pimpl->_charge; }

double ionic_group::ionion_energy_q() const
{
    //                                    q  // * number of protons
    if (_res.name() == "LYS") return 0.640; // * 81;
    if (_res.name() == "ASP") return 0.380; // * 95;
    if (_res.name() == "HIS") return 0.380; // * 83;
    if (_res.name() == "ARG") return 0.260; // * 77;
    if (_res.name() == "GLU") return 0.635; // * 69;

    // TODO exception it should not happen
    return 0;
}

string ionic_group::name() const
{
    return getNameFromAtoms(pimpl->_atoms);
}
