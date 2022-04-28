#include <string>
#include <stdexcept>

using namespace std;

#pragma region VDW

double opslVdw[][3] = {
        {0.000,  0.000, 0.000},//Index 1-based so [0] is a dummy
        {0.500,  3.750, 0.105},
        {-0.500, 2.960, 0.210},
        {-0.570, 3.250, 0.170},
        {0.370,  0.000, 0.000},
        {0.200,  3.800, 0.118},
        {0.200,  3.800, 0.080},
        {0.000,  3.910, 0.160},
        {0.000,  3.850, 0.080},
        {0.000,  3.905, 0.118},
        {0.000,  3.905, 0.175},
        {0.000,  3.750, 0.110},
        {-0.850, 3.250, 0.170},
        {0.425,  0.000, 0.000},
        {0.285,  3.800, 0.080},
        {0.285,  3.800, 0.118},
        {-0.100, 3.905, 0.118},
        {0.700,  3.750, 0.105},
        {-0.800, 2.960, 0.210},
        {0.310,  3.905, 0.118},
        {-0.300, 3.250, 0.170},
        {0.330,  0.000, 0.000},
        {0.265,  3.905, 0.118},
        {-0.700, 3.070, 0.170},
        {0.435,  0.000, 0.000},
        {0.265,  3.850, 0.080},
        {0.265,  3.750, 0.110},
        {0.310,  3.800, 0.118},
        {0.100,  3.800, 0.118},
        {0.310,  3.800, 0.080},
        {0.100,  3.800, 0.080},
        {0.180,  3.905, 0.118},
        {-0.450, 3.550, 0.250},
        {0.270,  0.000, 0.000},
        {0.235,  3.800, 0.118},
        {-0.470, 3.550, 0.250},
        {0.235,  3.800, 0.170},
        {0.300,  3.800, 0.118},
        {-0.300, 3.550, 0.250},
        {0.200,  3.800, 0.170},
        {-0.570, 3.250, 0.170},
        {0.420,  0.000, 0.000},
        {-0.490, 3.250, 0.170},
        {0.410,  3.750, 0.145},
        {0.100,  3.750, 0.145},
        {0.130,  3.750, 0.145},
        {-0.540, 3.250, 0.170},
        {0.460,  0.000, 0.000},
        {0.500,  3.750, 0.145},
        {0.330,  3.750, 0.145},
        {-0.055, 3.750, 0.145},
        {-0.800, 3.250, 0.170},
        {0.460,  0.000, 0.000},
        {0.640,  2.250, 0.050},
        {-0.700, 3.250, 0.170},
        {0.440,  0.000, 0.000},
        {0.310,  3.905, 0.118},
        {0.070,  3.905, 0.118},
        {0.550,  3.750, 0.105},
        {-0.450, 2.960, 0.210},
        {0.250,  3.800, 0.080},
        {0.250,  3.800, 0.118},
        {-0.400, 3.000, 0.170},
        {0.250,  3.800, 0.170},
        {0.200,  3.800, 0.050},
        {0.000,  3.960, 0.145},
};

int get_vdw_opsl_values_index(string const& residue, string const& atom, string const& element)
{
    if (residue == "GLY") {
        if (atom == "CA") return 5;
        if (element == "N") return 3;
        if (element == "C") return 1;
        if (element == "O") return 2;
    }
    if (residue == "PRO") {
        if (atom == "CA") return 14;
        if (atom == "CB") return 9;
        if (atom == "CG") return 9;
        if (atom == "CD") return 15;
        if (element == "N") return 3;
        if (element == "C") return 1;
        if (element == "O") return 2;
    }
    if (residue == "ALA") {
        if (atom == "CA") return 6;
        if (atom == "CB") return 7;
        if (element == "N") return 3;
        if (element == "C") return 1;
        if (element == "O") return 2;
    }
    if (residue == "AIB") {
        if (atom == "CA") return 64;
        if (atom == "CB") return 65;
        if (element == "N") return 3;
        if (element == "C") return 1;
        if (element == "O") return 2;
    }

    //First column
    if (residue == "ILE") {
        if (atom == "CB") return 8;
        if (atom == "CG") return 7;
        if (atom == "CD") return 10;
    }
    if (residue == "SER") {
        if (atom == "CB") return 22;
        if (atom == "OG") return 23;
    }
    if (residue == "THR") {
        if (atom == "CB") return 25;
        if (atom == "OG") return 23;
        if (atom == "CG") return 7;
    }
    if (residue == "TYR") {
        if (atom == "CB") return 9;
        if (atom == "CG") return 11;
        if (atom == "CD1") return 11;
        if (atom == "CD2") return 11;
        if (atom == "CE1") return 11;
        if (atom == "CE2") return 11;
        if (atom == "CZ") return 26;
        if (atom == "OH") return 23;
    }
    if (residue == "ASN") {
        if (atom == "CB") return 9;
        if (atom == "CG") return 1;
        if (atom == "OD1") return 2;
        if (atom == "ND2") return 12;
    }
    if (residue == "ASP") {
        if (atom == "CB") return 16;
        if (atom == "CG") return 17;
        if (atom == "OD1") return 18;
        if (atom == "OD2") return 18;
    }
    if (residue == "HIS") {
        if (atom == "CB") return 9;
        if (atom == "CG") return 45;
        if (atom == "ND1") return 40;
        if (atom == "CD2") return 44;
        if (atom == "CE1") return 43;
        if (atom == "NE2") return 42;
    }
    if (residue == "TRP") {
        if (atom == "CB") return 9;
        if (atom == "CG") return 50;
        if (atom == "CD1") return 45;
        if (atom == "CD2") return 50;
        if (atom == "NE1") return 40;
        if (atom == "CE2") return 45;
        if (atom == "CE3") return 11;
        if (atom == "CZ2") return 11;
        if (atom == "CZ3") return 11;
        if (atom == "CH2") return 11;
    }
    if (residue == "LYS") {
        if (atom == "CB") return 9;
        if (atom == "CG") return 9;
        if (atom == "CD") return 9;
        if (atom == "CE") return 19;
        if (atom == "NZ") return 20;
    }

    //Second column
    if (residue == "VAL") {
        if (atom == "CB") return 8;
        if (atom == "CG1") return 7;
        if (atom == "CG2") return 7;
    }
    if (residue == "LEU") {
        if (atom == "CB") return 9;
        if (atom == "CG") return 8;
        if (atom == "CD1") return 7;
        if (atom == "CD2") return 7;
    }
    if (residue == "PHE") {
        if (atom == "CB") return 9;
        if (atom == "CG") return 11;
        if (atom == "CD1") return 11;
        if (atom == "CD2") return 11;
        if (atom == "CE1") return 11;
        if (atom == "CE2") return 11;
        if (atom == "CZ") return 11;
    }
    if (residue == "CYS") {
        if (atom == "CB") return 9;
        if (atom == "SG") return 11;
    }
    if (residue == "MET") {
        if (atom == "CB") return 9;
        if (atom == "CG") return 34;
        if (atom == "SD") return 35;
        if (atom == "CE") return 36;
    }
    if (residue == "HIP") {
        if (atom == "CB") return 9;
        if (atom == "CG") return 49;
        if (atom == "ND1") return 46;
        if (atom == "CD2") return 49;
        if (atom == "CE1") return 48;
        if (atom == "NE2") return 46;
    }
    if (residue == "GLN") {
        if (atom == "CB") return 9;
        if (atom == "CG") return 9;
        if (atom == "CD") return 1;
        if (atom == "OE1") return 2;
        if (atom == "NE2") return 12;
        if (atom == "1HE2") return 13;
        if (atom == "2HE2") return 13;
    }
    if (residue == "GLU") {
        if (atom == "CB") return 9;
        if (atom == "CG") return 16;
        if (atom == "CD") return 17;
        if (atom == "OE1") return 18;
    }

    if (residue == "ARG") {
        if (atom == "CB") return 9;
        if (atom == "CG") return 57;
        if (atom == "CD") return 56;
        if (atom == "NE") return 54;
        if (atom == "CZ") return 53;
        if (atom == "NH1") return 51;
        if (atom == "NH2") return 51;
    }

    return -1; //Not present
}

bool has_vdw_opsl_values(string const& residue, string const& atom, string const& element)
{
    return get_vdw_opsl_values_index(residue, atom, element) > 0;
}


double* get_vdw_opsl_values(string const& residue, string const& atom, string const& element)
{ //Atom_name search must be before element search

    int index = get_vdw_opsl_values_index(residue, atom, element);
    if(index > 0)
        return opslVdw[index];
    else
        throw std::invalid_argument("get_vdw_opsl_values: residue " + residue + ", atom " + atom + ", element " + element + " unsupported");
}

#pragma endregion 