#include "ns_chemical_entity.h"

#include <memory>

#include "energy.h"

#include "private/impl_chemical_entity.h"

using std::vector, std::array, std::string, std::unique_ptr, std::make_unique, std::to_string, std::invalid_argument;

using chemical_entity::aminoacid, chemical_entity::atom, chemical_entity::ring, chemical_entity::ionic_group;

string join_strings(std::vector<std::string> const& values, string const& delimiter)
{
    string out;
    for (string const& value: values)
    {
        if (!out.empty())
            out += delimiter;
        out += value;
    }

    return out;
}

string get_name_from_atoms(vector<atom> const& atoms, string const& delimiter = ":")
{
    vector<string> atoms_name;
    for (auto const& atom: atoms)
        atoms_name.push_back(atom.get_name());

    sort(atoms_name.begin(), atoms_name.end());
    return join_strings(atoms_name, delimiter);
}

vector<atom> const& aminoacid::get_atoms() const
{ return _pimpl->atoms; }

string const& aminoacid::get_protein_name() const
{ return _pimpl->protein_name; }

std::optional<atom> const& aminoacid::get_alpha_carbon() const
{ return _pimpl->alpha_carbon; }

std::optional<atom> const& aminoacid::get_beta_carbon() const
{ return _pimpl->beta_carbon; }

std::optional<ring> const& aminoacid::get_primary_ring() const
{ return _pimpl->primary_ring; }

std::optional<ring> const& aminoacid::get_secondary_ring() const
{ return _pimpl->secondary_ring; }

std::optional<ionic_group> const& aminoacid::get_positive_ionic_group() const
{ return _pimpl->positive_ionic_group; }

std::optional<ionic_group> const& aminoacid::get_negative_ionic_group() const
{ return _pimpl->negative_ionic_group; }

string const& aminoacid::get_name() const
{ return _pimpl->name; }

string const& aminoacid::get_chain_id() const
{ return _pimpl->chain_id; }

string const& aminoacid::get_id() const
{ return _pimpl->id; }

int aminoacid::get_sequence_number() const
{ return _pimpl->sequence_number; }

bool aminoacid::operator==(aminoacid const& rhs) const
{ return _pimpl->id == rhs._pimpl->id; }

bool aminoacid::operator!=(aminoacid const& rhs) const
{ return !(*this == rhs); }

bool aminoacid::satisfies_minimum_sequence_separation(aminoacid const& other, int minimum_separation) const
{
    if (*this == other)
        return false;

    return get_chain_id() != other.get_chain_id() || abs(get_sequence_number() - other.get_sequence_number()) >= minimum_separation;
}

aminoacid::operator rin::node() const
{ return rin::node(*this); }

string aminoacid::get_secondary_structure_id() const
{ return _pimpl->secondary_structure_name; }

array<double, 3> center_of_mass(vector<atom> const& atoms)
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

// returning an optional of the error message should be way faster than try/catching
// exceptions (since a soon-to-be added option will give the possibility to ignore
// the exception and keep running)
std::optional<std::string> assert_atom_group_correctness(
    aminoacid const& residue, gemmi::Model const& model, vector<string> const& expected_atom_names, vector<atom> const& actual_atoms, string const& group_name)
{
    vector<string> actual_atom_names;
    actual_atom_names.reserve(actual_atoms.size());

    for(auto const& atom : actual_atoms)
        actual_atom_names.push_back(atom.get_name());

    for(auto const& expected_atom : expected_atom_names)
    {
        if (find(actual_atom_names.begin(), actual_atom_names.end(), expected_atom) == actual_atom_names.end())
        {
            auto const expected_atoms_str = join_strings(expected_atom_names, ",");
            auto const actual_atoms_str = get_name_from_atoms(actual_atoms, ",");

            string msg = "malformed " + group_name + " in model " + model.name;
            msg += " residue " + residue.get_id();
            msg += " expected atoms={";
            msg += expected_atoms_str;
            msg += "} actual atoms={";
            msg += actual_atoms_str;
            msg += "}";

            // an invalid atom group happened
            return msg;
        }
    }

    // no errors here
    return std::nullopt;
}

chemical_entity::aminoacid::aminoacid(
    gemmi::Residue const& residue,
    gemmi::Chain const& chain,
    gemmi::Model const& model,
    gemmi::Structure const& protein) :
    _pimpl(std::make_shared<impl>())
{
    _pimpl->name = residue.name;
    _pimpl->sequence_number = residue.seqid.num.value;
    _pimpl->chain_id = chain.name;

    _pimpl->id = _pimpl->chain_id + ":" + to_string(_pimpl->sequence_number) + ":_:" + _pimpl->name;

    // discover if this has 0, 1 or 2 aromatic rings
    int n_of_rings = 0;
    vector<string> ring1_names, ring2_names;

    if (_pimpl->name == "HIS")
    {
        // HIS has a 5-atoms ring
        ring1_names = {"CD2", "CE1", "CG", "ND1", "NE2"};
        n_of_rings = 1;
    }
    else if (_pimpl->name == "PHE" || _pimpl->name == "TYR")
    {
        // PHE, TYR have a 6-atoms ring
        ring1_names = {"CD1", "CD2", "CE1", "CE2", "CG", "CZ"};
        n_of_rings = 1;
    }
    else if (_pimpl->name == "TRP")
    {
        // TRP has both a 6-atoms ring and a 5-atoms ring
        ring1_names = {"CD2", "CE2", "CE3", "CH2", "CZ2", "CZ3"};
        ring2_names = {"CD2", "CE2", "CD1", "CG", "NE1"};

        n_of_rings = 2;
    }

    // discover ionic group of this aminoacid
    int charge = 0;
    vector<string> ionic_group_names;
    if (residue.name == "HIS")
    {
        charge = 1;
        ionic_group_names = {"CG", "CD2", "CE1", "ND1", "NE2"};
    }
    else if (residue.name == "ARG")
    {
        charge = 1;
        ionic_group_names = {"CZ", "NH2", "NH1", "NE"};
    }
    else if (residue.name == "LYS")
    {
        charge = 1;
        ionic_group_names = {"NZ"};
    }
    else if (residue.name == "GLU")
    {
        charge = -1;
        ionic_group_names = {"CD", "OE1", "OE2"};
    }
    else if (residue.name == "ASP")
    {
        charge = -1;
        ionic_group_names = {"CG", "OD1", "OD2"};
    }
    vector<atom> ring1, ring2, ionic_group_atoms;

    for (auto const& record: residue.atoms)
    {
        auto atom = chemical_entity::atom{record, *this};
        _pimpl->atoms.push_back(atom);

        if (atom.get_name() == "CA")
            _pimpl->alpha_carbon = atom;

        else if (atom.get_name() == "CB")
            _pimpl->beta_carbon = atom;

        if (n_of_rings >= 1 && find(ring1_names.begin(), ring1_names.end(), atom.get_name()) != ring1_names.end())
            ring1.push_back(atom);

        if (n_of_rings == 2 && find(ring2_names.begin(), ring2_names.end(), atom.get_name()) != ring2_names.end())
            ring2.push_back(atom);

        if (charge != 0 && find(ionic_group_names.begin(), ionic_group_names.end(), atom.get_name()) != ionic_group_names.end())
            ionic_group_atoms.push_back(atom);
    }

    _pimpl->position = center_of_mass(get_atoms());

    if (n_of_rings >= 1)
    {
        //assert_atom_group_correctness(*this, model, ring1_names, ring1, "ring");
        _pimpl->primary_ring = ring(ring1, *this);
    }
    if (n_of_rings == 2)
    {
       //assert_atom_group_correctness(*this, model, ring2_names, ring2, "ring");
       _pimpl->secondary_ring = ring(ring2, *this);
    }

    if (charge == 1)
    {
        //assert_atom_group_correctness(*this, model, ionic_group_names, ionic_group_atoms, "ionic group");
        _pimpl->positive_ionic_group = ionic_group(ionic_group_atoms, 1, *this);
    }
    else if (charge == -1)
    {
        //assert_atom_group_correctness(*this, model, ionic_group_names, ionic_group_atoms, "ionic group");
        _pimpl->negative_ionic_group = ionic_group(ionic_group_atoms, -1, *this);
    }

    _pimpl->protein_name = protein.name;
    _pimpl->secondary_structure_name = cfg::graphml::none;
}

aminoacid::aminoacid(
    gemmi::Residue const& residue,
    gemmi::Chain const& chain,
    gemmi::Model const& model,
    gemmi::Structure const& protein,
    std::optional<std::variant<gemmi::Helix, gemmi::Sheet::Strand>> const& secondary_structure) :
    aminoacid(residue, chain, model, protein)
{
    if (!secondary_structure.has_value())
        _pimpl->secondary_structure_name = "LOOP";
    else if (std::holds_alternative<gemmi::Helix>(*secondary_structure))
        _pimpl->secondary_structure_name = "HELIX";
    else if (std::holds_alternative<gemmi::Sheet::Strand>(*secondary_structure))
        _pimpl->secondary_structure_name = "SHEET";
}

aminoacid aminoacid::component::get_residue() const
{
    // information-less aminoacid
    aminoacid res;

    // restore all of its information
    res._pimpl = _res_pimpl.lock();
    return res;
}

aminoacid::aminoacid() = default;

aminoacid::~aminoacid() = default;

std::array<double, 3> const& chemical_entity::aminoacid::get_position() const
{ return _pimpl->position; }

atom::atom(gemmi::Atom const& record, aminoacid const& res) :
    kdpoint<3>({record.pos.x, record.pos.y, record.pos.z}),
    component(res),
    _pimpl{new impl{record}}
{}

atom::~atom() = default;

string const& atom::get_name() const
{ return _pimpl->record.name; }

string atom::get_symbol() const
{ return gemmi::element_uppercase_name(_pimpl->record.element.elem); }

double atom::get_temp_factor() const
{ return _pimpl->record.b_iso; }

int atom::get_charge() const
{
    auto c = _pimpl->record.charge;
    return c > 0 ? 1 : c < 0 ? -1 : 0;
}

bool atom::is_hydrogen() const
{ return _pimpl->record.is_hydrogen(); }

bool atom::is_main_chain() const
{
    auto names = {"C", "O", "H", "HA", "N"};
    return find(names.begin(), names.end(), get_name()) != names.end();
}

int atom::get_atom_number() const
{ return _pimpl->record.serial; }

double atom::get_mass() const
{
    auto element = get_symbol();
    if (element == "H") return 1.008;
    if (element == "C") return 12.011;
    if (element == "N") return 14.007;
    if (element == "O") return 15.994;
    if (element == "S") return 32.065;

    return _pimpl->record.element.weight();
    //throw std::invalid_argument("atom::get_mass(): unsupported element " + element);
}

double atom::get_vdw_radius() const
{
    string element = get_symbol();
    if (element == "S") return 1.89;
    if (element == "C") return 1.77;
    if (element == "O") return 1.55;
    if (element == "N") return 1.60;

    return _pimpl->record.element.vdw_r();
    //throw std::invalid_argument("atom::get_vdw_radius(): unsupported element " + element);
}

bool atom::is_cation() const
{
    auto res_name = get_residue().get_name();

    return (res_name == "LYS" && get_name() == "NZ")
        || (res_name == "ARG" && get_name() == "NH2")
        || (res_name == "HIS" && get_name() == "ND1");
}

bool atom::in_positive_ionic_group() const
{
    auto res_name = get_residue().get_name();
    auto name = get_name();

    if (res_name == "HIS")
        return name == "CG" || name == "CD2" || name == "CE1" || name == "ND1" || name == "NE2";

    else if (res_name == "ARG")
        return name == "CZ" || name == "NH2" || name == "NH1" || name == "NE";

    else if (res_name == "LYS")
        return name == "NZ";

    return false;
}

bool chemical_entity::atom::in_negative_ionic_group() const
{
    auto res_name = get_residue().get_name();
    auto name = get_name();

    if (res_name == "GLU")
        return name == "CD" || name == "OE1" || name == "OE2";

    else if (res_name == "ASP")
        return name == "CG" || name == "OD1" || name == "OD2";

    return false;
}

bool atom::is_hydrogen_donor() const
{
    auto res_name = get_residue().get_name();
    auto name = get_name();
    return (res_name == "ARG" && (name == "NH1" || name == "NH2" || name == "NE"))
        || (res_name == "ASN" && name == "ND2")
        || (res_name == "GLN" && name == "NE2")
        || (res_name == "HIS" && (name == "NE2" || name == "ND1"))
        || (res_name == "LYS" && name == "NZ")
        || (res_name == "SER" && name == "OG")
        || (res_name == "THR" && name == "OG1")
        || (res_name == "TRP" && name == "NE1")
        || (res_name == "TYR" && name == "OH")
        || name == "NH"
        || name == "N";
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
    auto name = get_name();

    return (res_name == "ASN" && name == "OD1")
        || (res_name == "ASP" && (name == "OD1" || name == "OD2"))
        || (res_name == "GLN" && name == "OE1")
        || (res_name == "GLU" && (name == "OE1" || name == "OE2"))
        || (res_name == "HIS" && (name == "ND1" || name == "NE2"))
        || (res_name == "SER" && name == "OG")
        || (res_name == "THR" && name == "OG1")
        || (res_name == "TYR" && name == "OH")
        //|| n == "C"
        || name == "O";
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
    auto name = get_name();
    auto element = get_symbol();

    return has_vdw_opsl_values(res_name, name, element);

    //return (res_name == "GLN" && (name == "NE1" || name == "OE1"))
    //    || (res_name == "ASN" && (name == "ND2" || name == "OD1"))
    //    || element == "C"
    //    || element == "S";
}

vector<atom> atom::get_attached_hydrogens() const
{
    vector<atom> hydrogens;
    auto const hydrogen_name_pattern = "H" + get_name().substr(1, get_name().size() - 1);
    for (auto const& atom : get_residue().get_atoms())
    {
        if (atom.is_hydrogen() && atom.get_name().find(hydrogen_name_pattern) != std::string::npos)
            hydrogens.push_back(atom);
    }

    return hydrogens;
}

ring::ring(vector<atom> const& atoms, aminoacid const& res) : kdpoint<3>({0, 0, 0}), component(res)
{
    auto tmp_pimpl = std::make_shared<impl>();

    // already taken care of in assert_ring_correctness
    //if (atoms.size() < 3) throw invalid_argument("rings should have at least 3 atoms");

    tmp_pimpl->atoms = atoms;

    _position = center_of_mass(atoms);

    // kudos to Giulio Marcolin for the following shortcut
    // it only deviates from a SVD best-fit method no more than 1-2°, on average
    array<double, 3> const v = (array<double, 3>) (atoms[0] - atoms[1]);
    array<double, 3> const w = (array<double, 3>) (atoms[2] - atoms[1]);
    tmp_pimpl->normal = geom::cross(v, w);

    _pimpl = tmp_pimpl;
}

ring::~ring() = default;

array<double, 3> const& ring::get_normal() const
{ return _pimpl->normal; }

bool ring::is_pication_candidate() const
{
    string name = get_residue().get_name();
    return name == "PHE" || name == "TYR" || (name == "TRP" && _pimpl->atoms.size() == 6);
}

double ring::get_distance_between_closest_atoms(ring const& other) const
{
    double minimum = _pimpl->atoms[0].distance(other._pimpl->atoms[0]);
    for (auto const& atom_1: _pimpl->atoms)
    {
        for (auto const& atom_2: other._pimpl->atoms)
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
{ return geom::d_angle<3>(_pimpl->normal, other._pimpl->normal); }

double ring::get_angle_between_normal_and_centers_joining(ring const& other) const
{
    std::array<double, 3> const centres_joining((std::array<double, 3>) (*this - other));
    return geom::d_angle<3>(_pimpl->normal, centres_joining);
}

string ring::get_name() const
{ return get_name_from_atoms(_pimpl->atoms); }

ionic_group::ionic_group(vector<atom> const& atoms, int const& charge, aminoacid const& res) :
    kdpoint<3>({0, 0, 0}), component(res), _pimpl{new impl{atoms, charge}}
{ _position = center_of_mass(atoms); }

ionic_group::~ionic_group() = default;

int ionic_group::get_charge() const
{ return _pimpl->charge; }

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
{ return get_name_from_atoms(_pimpl->atoms); }
