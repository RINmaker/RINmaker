#pragma once

#include <string>

/// <summary>
/// It return the value needed to calculate the vdw energy
/// </summary>
/// <param name="residue">Aminoacid three-letters code</param>
/// <param name="atom_name">PDB atom name</param>
/// <param name="element_name">Element name</param>
/// <returns>{q, sigma, epsilon (kcal/mol)}</returns>
double *get_vdw_opsl_values(std::string residue, std::string atom_name, std::string element_name);