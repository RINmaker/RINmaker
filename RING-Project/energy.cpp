#include <string>
using namespace std;

#pragma region VDW

double opslVdw[][3] = {
{ 0.000, 0.000, 0.000},//Index 1-based so [0] is a dummy
{ 0.500, 3.750, 0.105},
{-0.500, 2.960, 0.210},
{-0.570, 3.250, 0.170},
{ 0.370, 0.000, 0.000},
{ 0.200, 3.800, 0.118},
{ 0.200, 3.800, 0.080},
{ 0.000, 3.910, 0.160},
{ 0.000, 3.850, 0.080},
{ 0.000, 3.905, 0.118},
{ 0.000, 3.905, 0.175},
{ 0.000, 3.750, 0.110},
{-0.850, 3.250, 0.170},
{ 0.425, 0.000, 0.000},
{ 0.285, 3.800, 0.080},
{ 0.285, 3.800, 0.118},
{-0.100, 3.905, 0.118},
{ 0.700, 3.750, 0.105},
{-0.800, 2.960, 0.210},
{ 0.310, 3.905, 0.118},
{-0.300, 3.250, 0.170},
{ 0.330, 0.000, 0.000},
{ 0.265, 3.905, 0.118},
{-0.700, 3.070, 0.170},
{ 0.435, 0.000, 0.000},
{ 0.265, 3.850, 0.080},
{ 0.265, 3.750, 0.110},
{ 0.310, 3.800, 0.118},
{ 0.100, 3.800, 0.118},
{ 0.310, 3.800, 0.080},
{ 0.100, 3.800, 0.080},
{ 0.180, 3.905, 0.118},
{-0.450, 3.550, 0.250},
{ 0.270, 0.000, 0.000},
{ 0.235, 3.800, 0.118},
{-0.470, 3.550, 0.250},
{ 0.235, 3.800, 0.170},
{ 0.300, 3.800, 0.118},
{-0.300, 3.550, 0.250},
{ 0.200, 3.800, 0.170},
{-0.570, 3.250, 0.170},
{ 0.420, 0.000, 0.000},
{-0.490, 3.250, 0.170},
{ 0.410, 3.750, 0.145},
{ 0.100, 3.750, 0.145},
{ 0.130, 3.750, 0.145},
{-0.540, 3.250, 0.170},
{ 0.460, 0.000, 0.000},
{ 0.500, 3.750, 0.145},
{ 0.330, 3.750, 0.145},
{-0.055, 3.750, 0.145},
{-0.800, 3.250, 0.170},
{ 0.460, 0.000, 0.000},
{ 0.640, 2.250, 0.050},
{-0.700, 3.250, 0.170},
{ 0.440, 0.000, 0.000},
{ 0.310, 3.905, 0.118},
{ 0.070, 3.905, 0.118},
{ 0.550, 3.750, 0.105},
{-0.450, 2.960, 0.210},
{ 0.250, 3.800, 0.080},
{ 0.250, 3.800, 0.118},
{-0.400, 3.000, 0.170},
{ 0.250, 3.800, 0.170},
{ 0.200, 3.800, 0.050},
{ 0.000, 3.960, 0.145},
};

double* get_vdw_opsl_values(string residue, string atom_name, string element_name)
{ //Atom_name search must be before element_name search

	auto aux = [residue, atom_name, element_name]() {
		if (residue == "GLY")
		{
			if (atom_name == "CA") return 5;
			if (element_name == "N") return 3;
			if (element_name == "C") return 1;
			if (element_name == "O") return 2;
		}
		if (residue == "PRO")
		{
			if (atom_name == "CA") return 14;
			if (atom_name == "CB") return 9;
			if (atom_name == "CG") return 9;
			if (atom_name == "CD") return 15;
			if (element_name == "N") return 3;
			if (element_name == "C") return 1;
			if (element_name == "O") return 2;
		}
		if (residue == "ALA")
		{
			if (atom_name == "CA") return 6;
			if (atom_name == "CB") return 7;
			if (element_name == "N") return 3;
			if (element_name == "C") return 1;
			if (element_name == "O") return 2;
		}
		if (residue == "AIB")
		{
			if (atom_name == "CA") return 64;
			if (atom_name == "CB") return 65;
			if (element_name == "N") return 3;
			if (element_name == "C") return 1;
			if (element_name == "O") return 2;
		}

		//First column
		if (residue == "ILE")
		{
			if (atom_name == "CB") return 8;
			if (atom_name == "CG") return 7;
			if (atom_name == "CD") return 10;
		}
		if (residue == "SER")
		{
			if (atom_name == "CB") return 22;
			if (atom_name == "OG") return 23;
		}
		if (residue == "THR")
		{
			if (atom_name == "CB") return 25;
			if (atom_name == "OG") return 23;
			if (atom_name == "CG") return 7;
		}
		if (residue == "TYR")
		{
			if (atom_name == "CB") return 9;
			if (atom_name == "CG") return 11;
			if (atom_name == "CD1") return 11;
			if (atom_name == "CD2") return 11;
			if (atom_name == "CE1") return 11;
			if (atom_name == "CE2") return 11;
			if (atom_name == "CZ") return 26;
			if (atom_name == "OH") return 23;
		}
		if (residue == "ASN")
		{
			if (atom_name == "CB") return 9;
			if (atom_name == "CG") return 1;
			if (atom_name == "OD1") return 2;
			if (atom_name == "ND2") return 12;
		}
		if (residue == "ASP")
		{
			if (atom_name == "CB") return 16;
			if (atom_name == "CG") return 17;
			if (atom_name == "OD1") return 18;
			if (atom_name == "OD2") return 18;
		}
		if (residue == "HIS")
		{
			if (atom_name == "CB") return 9;
			if (atom_name == "CG") return 45;
			if (atom_name == "ND1") return 40;
			if (atom_name == "CD2") return 44;
			if (atom_name == "CE1") return 43;
			if (atom_name == "NE2") return 42;
		}
		if (residue == "TRP")
		{
			if (atom_name == "CB") return 9;
			if (atom_name == "CG") return 50;
			if (atom_name == "CD1") return 45;
			if (atom_name == "CD2") return 50;
			if (atom_name == "NE1") return 40;
			if (atom_name == "CE2") return 45;
			if (atom_name == "CE3") return 11;
			if (atom_name == "CZ2") return 11;
			if (atom_name == "CZ3") return 11;
			if (atom_name == "CH2") return 11;
		}
		if (residue == "LYS")
		{
			if (atom_name == "CB") return 9;
			if (atom_name == "CG") return 9;
			if (atom_name == "CD") return 9;
			if (atom_name == "CE") return 19;
			if (atom_name == "NZ") return 20;
		}

		//Second column
		if (residue == "VAL")
		{
			if (atom_name == "CB") return 8;
			if (atom_name == "CG1") return 7;
			if (atom_name == "CG2") return 7;
		}
		if (residue == "LEU")
		{
			if (atom_name == "CB") return 9;
			if (atom_name == "CG") return 8;
			if (atom_name == "CD1") return 7;
			if (atom_name == "CD2") return 7;
		}
		if (residue == "PHE")
		{
			if (atom_name == "CB") return 9;
			if (atom_name == "CG") return 11;
			if (atom_name == "CD1") return 11;
			if (atom_name == "CD2") return 11;
			if (atom_name == "CE1") return 11;
			if (atom_name == "CE2") return 11;
			if (atom_name == "CZ") return 11;
		}
		if (residue == "CYS")
		{
			if (atom_name == "CB") return 9;
			if (atom_name == "SG") return 11;
		}
		if (residue == "MET")
		{
			if (atom_name == "CB") return 9;
			if (atom_name == "CG") return 34;
			if (atom_name == "SD") return 35;
			if (atom_name == "CE") return 36;
		}
		if (residue == "HIP")
		{
			if (atom_name == "CB") return 9;
			if (atom_name == "CG") return 49;
			if (atom_name == "ND1") return 46;
			if (atom_name == "CD2") return 49;
			if (atom_name == "CE1") return 48;
			if (atom_name == "NE2") return 46;
		}
		if (residue == "GLN")
		{
			if (atom_name == "CB") return 9;
			if (atom_name == "CG") return 9;
			if (atom_name == "CD") return 1;
			if (atom_name == "OE1") return 2;
			if (atom_name == "NE2") return 12;
			if (atom_name == "1HE2") return 13;
			if (atom_name == "2HE2") return 13;
		}
		if (residue == "GLU")
		{
			if (atom_name == "CB") return 9;
			if (atom_name == "CG") return 16;
			if (atom_name == "CD") return 17;
			if (atom_name == "OE1") return 18;
		}

		if (residue == "ARG")
		{
			if (atom_name == "CB") return 9;
			if (atom_name == "CG") return 57;
			if (atom_name == "CD") return 56;
			if (atom_name == "NE") return 54;
			if (atom_name == "CZ") return 53;
			if (atom_name == "NH1") return 51;
			if (atom_name == "NH2") return 51;
		}

		return 0;
	};

	return opslVdw[aux()];
}

#pragma endregion 