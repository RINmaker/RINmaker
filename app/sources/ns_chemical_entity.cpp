#include "ns_chemical_entity.h"

#include <memory>

#include "energy.h"

#include "prelude.h"
#include "utils.h"

#include "private/impl_chemical_entity.h"

using std::vector, std::array, std::string, std::unique_ptr, std::make_unique, std::to_string, std::invalid_argument;

using chemical_entity::aminoacid, chemical_entity::atom, chemical_entity::ring, chemical_entity::ionic_group;

string getNameFromAtoms(vector<atom> const& atoms, string const& delimiter = ":")
{
    vector<string> atoms_name;
    for (auto const& a: atoms)
        atoms_name.push_back(a.get_name());

    sort(atoms_name.begin(), atoms_name.end());

    return joinStrings(atoms_name, delimiter);
}

vector<atom> const& aminoacid::get_atoms() const
{
    return pimpl->_atoms;
}

string const& aminoacid::get_protein_name() const
{ return pimpl->pdb_name; }

std::optional<atom> const& aminoacid::get_alpha_carbon() const
{ return pimpl->alpha_carbon; }

std::optional<atom> const& aminoacid::get_beta_carbon() const
{ return pimpl->beta_carbon; }

std::optional<ring> const& aminoacid::get_primary_ring() const
{ return pimpl->primary_ring; }

std::optional<ring> const& aminoacid::get_secondary_ring() const
{ return pimpl->secondary_ring; }

std::optional<ionic_group> const& aminoacid::get_positive_ionic_group() const
{ return pimpl->positive_ionic_group; }

std::optional<ionic_group> const& aminoacid::get_negative_ionic_group() const
{ return pimpl->negative_ionic_group; }

string const& aminoacid::get_name() const
{ return pimpl->name; }

string const& aminoacid::get_chain_id() const
{ return pimpl->chain_id; }

string const& aminoacid::get_id() const
{ return pimpl->id; }

int aminoacid::get_sequence_number() const
{ return pimpl->sequence_number; }

bool aminoacid::operator==(aminoacid const& rhs) const
{ return pimpl->id == rhs.pimpl->id; }

bool aminoacid::operator!=(aminoacid const& rhs) const
{ return !(*this == rhs); }

bool aminoacid::satisfies_minimum_sequence_separation(aminoacid const& other, int minimum_separation) const
{
    if (*this == other)
    {
        return false;
    }

    if (get_chain_id() != other.get_chain_id())
    {
        return true;
    }

    return abs(get_sequence_number() - other.get_sequence_number()) >= minimum_separation;
}

aminoacid::operator rin::node() const
{ return rin::node(*this); }

string aminoacid::get_secondary_structure_id() const
{
    return pimpl->secondary_structure_name;
}

array<double, 3> centre_of_mass(vector<atom> const& atoms)
{
    double mass{0.0};
    array<double, 3> centroid{0.0, 0.0, 0.0};
    for (auto const& a: atoms)
    {
        double const m = a.get_mass();
        for (size_t i = 0; i < 3; ++i)
            centroid[i] += a[i] * m;
        mass += m;
    }
    for (size_t i = 0; i < 3; ++i)
        centroid[i] /= mass;

    return centroid;
}

void assert_ring_correctness(
        string const& name, uint32_t line_number, vector<string> const& expected_atoms, vector<atom> const& found_atoms)
{
    if (expected_atoms.size() != found_atoms.size())
    {
        string expected_atoms_str = joinStrings(expected_atoms, ",");
        string found_atoms_str = getNameFromAtoms(found_atoms, ",");
        string exception_description =
                "line number: " + std::to_string(line_number) + ", aminoacid: " + name +
                " - expected aromatic ring: " + expected_atoms_str + ", found: " + found_atoms_str;
        throw std::invalid_argument(exception_description);
    }
};

chemical_entity::aminoacid::aminoacid(
    gemmi::Residue const& residue,
    gemmi::Chain const& chain,
    gemmi::Model const& model,
    gemmi::Structure const& protein) :
    pimpl(std::make_shared<impl>())
{
    pimpl->name = residue.name;
    pimpl->sequence_number = residue.seqid.num.value;
    pimpl->chain_id = chain.name;

    pimpl->id = pimpl->chain_id + ":" + to_string(pimpl->sequence_number) + ":_:" + pimpl->name;

    // discover if this has 0, 1 or 2 aromatic rings
    int n_of_rings = 0;
    vector<string> patterns_1, patterns_2;

    if (pimpl->name == "HIS")
    {
        // HIS has a 5-atoms ring
        patterns_1 = {"CD2", "CE1", "CG", "ND1", "NE2"};
        n_of_rings = 1;
    }
    else if (pimpl->name == "PHE" || pimpl->name == "TYR")
    {
        // PHE, TYR have a 6-atoms ring
        patterns_1 = {"CD1", "CD2", "CE1", "CE2", "CG", "CZ"};
        n_of_rings = 1;
    }
    else if (pimpl->name == "TRP")
    {
        // TRP has both a 6-atoms ring and a 5-atoms ring
        patterns_1 = {"CD2", "CE2", "CE3", "CH2", "CZ2", "CZ3"};
        patterns_2 = {"CD2", "CE2", "CD1", "CG", "NE1"};

        n_of_rings = 2;
    }

    // note: are they mutually exclusive? should be addressed
    vector<atom> positive, negative;
    vector<atom> ring_1, ring_2;

    for (auto const& record: residue.atoms)
    {
        auto atom = chemical_entity::atom{record, *this};
        pimpl->_atoms.push_back(atom);

        if (atom.get_name() == "CA")
            pimpl->alpha_carbon = atom;

        else if (atom.get_name() == "CB")
            pimpl->beta_carbon = atom;

        if (n_of_rings >= 1 && find(patterns_1.begin(), patterns_1.end(), atom.get_name()) != patterns_1.end())
            ring_1.push_back(atom);

        if (n_of_rings == 2 && find(patterns_2.begin(), patterns_2.end(), atom.get_name()) != patterns_2.end())
            ring_2.push_back(atom);

        if (atom.in_positive_ionic_group())
            positive.push_back(atom);

        else if (atom.in_negative_ionic_group())
            negative.push_back(atom);
    }

    pimpl->pos = centre_of_mass(get_atoms());

    if (n_of_rings >= 1)
    {
        // fixme find a way to locate line numbers from gemmi
        assert_ring_correctness(pimpl->name, -1, patterns_1, ring_1);
        pimpl->primary_ring = ring(ring_1, *this);
    }
    if (n_of_rings == 2)
    {
        // fixme find a way to locate line numbers from gemmi
        assert_ring_correctness(pimpl->name, -1, patterns_2, ring_2);
        pimpl->secondary_ring = ring(ring_2, *this);
    }

    if (!positive.empty())
        pimpl->positive_ionic_group = ionic_group(positive, 1, *this);

    if (!negative.empty())
        pimpl->negative_ionic_group = ionic_group(negative, -1, *this);

    pimpl->pdb_name = protein.name;

    if (protein.helices.empty() && protein.sheets.empty())
        pimpl->secondary_structure_name = "NONE";
    else
    {
        pimpl->secondary_structure_name = "LOOP";

        // todo retrieve secondary structure
        //if (model.find_cra(protein.helices.front().res).chain->find_residue(residue) != nullptr
        //||  model.find_cra(protein.helices.front().start).chain->find_residue(residue);
    }
}

aminoacid aminoacid::component::get_residue() const
{
    // information-less aminoacid
    aminoacid res;

    // restore all of its information
    // TODO check and throw exception?
    // to me it is redundant, because weak_ptr already throws iff not valid when calling .lock()
    res.pimpl = res_impl.lock();
    return res;
}


aminoacid::aminoacid() = default;

aminoacid::~aminoacid() = default;

std::array<double, 3> const& chemical_entity::aminoacid::get_position() const
{
    return pimpl->pos;
}

atom::atom(gemmi::Atom const& record, aminoacid const& res) :
    kdpoint<3>({record.pos.x, record.pos.y, record.pos.z}),
    component(res),
    pimpl{new impl{record}}
{}

atom::~atom() = default;

string const& atom::get_name() const
{ return pimpl->record.name; }

string atom::get_symbol() const
{ return gemmi::element_uppercase_name(pimpl->record.element.elem); }

double atom::get_temp_factor() const
{ return pimpl->record.b_iso; }

int atom::get_charge() const
{
    auto c = pimpl->record.charge;
    return c > 0 ? 1 : c < 0 ? -1 : 0;
}

bool atom::is_hydrogen() const
{ return pimpl->record.is_hydrogen(); }

bool atom::is_main_chain() const
{
    return
        this->get_name() == "C" ||
                this->get_name() == "O" ||
                this->get_name() == "H" ||
                this->get_name() == "HA" ||
                this->get_name() == "N";
}


int atom::get_atom_number() const
{ return pimpl->record.serial; }

double atom::get_mass() const
{
    auto element = get_symbol();
    if (element == "H") return 1.008;
    if (element == "C") return 12.011;
    if (element == "N") return 14.007;
    if (element == "O") return 15.994;
    if (element == "S") return 32.065;

    throw std::invalid_argument("atom::get_mass(): unsupported element " + element);
}

double atom::get_vdw_radius() const
{
    string element = get_symbol();
    if (element == "S") return 1.89;
    if (element == "C") return 1.77;
    if (element == "O") return 1.55;
    if (element == "N") return 1.60;

    throw std::invalid_argument("atom::get_vdw_radius(): unsupported element " + element);
}

bool atom::is_cation() const
{
    auto res_name = get_residue().get_name();

    return (res_name == "LYS" && get_name() == "NZ") ||
           (res_name == "ARG" && get_name() == "NH2") ||
           (res_name == "HIS" && get_name() == "ND1");
}

bool atom::in_positive_ionic_group() const
{
    auto res_name = get_residue().get_name();

    if (res_name == "HIS")
    {
        return get_name() == "CG" || get_name() == "CD2" || get_name() == "CE1" || get_name() == "ND1" || get_name() == "NE2";
    }
    else if (res_name == "ARG")
    {
        return get_name() == "CZ" || get_name() == "NH2";
    }
    else if (res_name == "LYS")
    {
        return get_name() == "NZ";
    }
    return false;
}

bool chemical_entity::atom::in_negative_ionic_group() const
{
    auto res_name = get_residue().get_name();
    auto n = get_name();

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

bool atom::is_hydrogen_donor() const
{
    auto res_name = get_residue().get_name();
    auto n = get_name();
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
            n == "NH" || n == "N";
}

int atom::how_many_hydrogen_can_donate() const
{
    if (is_hydrogen_donor())
    {
        std::string res_name = get_residue().get_name();
        std::string n = get_name();
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

bool atom::is_hydrogen_acceptor() const
{
    auto res_name = get_residue().get_name();
    auto n = get_name();

    return
            (res_name == "ASN" && n == "OD1") ||
            (res_name == "ASP" && (n == "OD1" || n == "OD2")) ||
            (res_name == "GLN" && n == "OE1") ||
            (res_name == "GLU" && (n == "OE1" || n == "OE2")) ||
            (res_name == "HIS" && (n == "ND1" || n == "NE2")) ||
            (res_name == "SER" && n == "OG") ||
            (res_name == "THR" && n == "OG1") ||
            (res_name == "TYR" && n == "OH") ||
            //n == "C" ||
            n == "O";
}

int atom::how_many_hydrogen_can_accept() const
{
    if (is_hydrogen_acceptor())
    {
        std::string res_name = get_residue().get_name();
        std::string n = get_name();
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

bool atom::is_vdw_candidate() const
{
    auto res_name = get_residue().get_name();
    auto n = get_name();
    auto en = get_symbol();

    return has_vdw_opsl_values(res_name, n, en);

    /*return
            (res_name == "GLN" && (n == "NE1" || n == "OE1")) ||
            (res_name == "ASN" && (n == "ND2" || n == "OD1")) ||
            en == "C" || en == "S";*/
}

vector<atom> atom::get_attached_hydrogens() const
{
    vector<atom> hydrogens;
    auto const hydrogen_name_pattern = "H" + get_name().substr(1, get_name().size() - 1);
    for (auto const& a : get_residue().get_atoms())
    {
        if (a.is_hydrogen() && prelude::match(a.get_name(), hydrogen_name_pattern))
            hydrogens.push_back(a);
    }

    return hydrogens;
}

ring::ring(vector<atom> const& atoms, aminoacid const& res) :
        kdpoint<3>({0, 0, 0}), component(res)
{
    auto tmp_pimpl = std::make_shared<impl>();

    if (atoms.size() < 3)
        throw invalid_argument("rings should have at least 3 atoms");

    tmp_pimpl->atoms = atoms;

    double sum_radii = 0;
    for (auto const& a: atoms)
        sum_radii += distance(a);
    tmp_pimpl->mean_radius = sum_radii / (double) atoms.size();

    _position = centre_of_mass(atoms);

    // kudos to Giulio Marcolin for the following shortcut
    // it only deviates from a SVD best-fit method no more than 1-2°, on average
    array<double, 3> const v = (array<double, 3>) (atoms[0] - atoms[1]);
    array<double, 3> const w = (array<double, 3>) (atoms[2] - atoms[1]);
    tmp_pimpl->normal = geom::cross(v, w);

    pimpl = tmp_pimpl;
}

ring::~ring() = default;

array<double, 3> const& ring::get_normal() const
{ return pimpl->normal; }

bool ring::is_pication_candidate() const
{
    string name = get_residue().get_name();
    return name == "PHE" || name == "TYR" || (name == "TRP" && pimpl->atoms.size() == 6);
}

double ring::get_distance_between_closest_atoms(ring const& other) const
{
    double minimum = pimpl->atoms[0].distance(other.pimpl->atoms[0]);
    for (auto const& atom_1: pimpl->atoms)
    {
        for (auto const& atom_2: other.pimpl->atoms)
        {
            double current = atom_1.distance(atom_2);
            if (current < minimum)
            {
                minimum = current;
            }
        }
    }

    return minimum;
}

double ring::get_angle_between_normals(ring const& other) const
{ return geom::d_angle<3>(pimpl->normal, other.pimpl->normal); }

double ring::get_angle_between_normal_and_centers_joining(ring const& other) const
{
    std::array<double, 3> const centres_joining((std::array<double, 3>) (*this - other));
    return geom::d_angle<3>(pimpl->normal, centres_joining);
}

string ring::get_name() const
{
    return getNameFromAtoms(pimpl->atoms);
}

ionic_group::ionic_group(vector<atom> const& atoms, int const& charge, aminoacid const& res) :
        kdpoint<3>({0, 0, 0}), component(res), pimpl{new impl{atoms, charge}}
{ _position = centre_of_mass(atoms); }

ionic_group::~ionic_group() = default;

int ionic_group::get_charge() const
{ return pimpl->charge; }

double ionic_group::get_ionion_energy_q() const
{
    string res_name = get_residue().get_name();
    //                                q  // * number of protons
    if (res_name == "LYS") return 0.640; // * 81;
    if (res_name == "ASP") return 0.380; // * 95;
    if (res_name == "HIS") return 0.380; // * 83;
    if (res_name == "ARG") return 0.260; // * 77;
    if (res_name == "GLU") return 0.635; // * 69;

    throw std::invalid_argument("ionic_group::get_ionion_energy_q(): unsupported residue " + res_name);
}

string ionic_group::get_name() const
{
    return getNameFromAtoms(pimpl->atoms);
}
