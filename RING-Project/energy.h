#pragma once

#include <string>
using namespace std;

/// <summary>
/// It return the value needed to calculate the vdw energy
/// </summary>
/// <param name="residue">Aminoacid three-letters code</param>
/// <param name="atom_name">PDB atom name</param>
/// <param name="element_name">Element name</param>
/// <returns>{q, sigma, epsilon (kcal/mol)}</returns>
double* get_vdw_opsl_values(string residue, string atom_name, string element_name);