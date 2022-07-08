#pragma once

#include <string>

/// It return the value needed to calculate the vdw energy
/// \param residue Aminoacid three-letters code
/// \param atom PDB atom name
/// \param element Element name
/// \return {q, sigma, epsilon (kcal/mol)}
double* get_vdw_opsl_values(std::string const& residue, std::string const& atom, std::string const& element);

/// Indicate if there is information to calculate the vdw energy
/// \param residue Aminoacid three-letters code
/// \param atom PDB atom name
/// \param element Element name
bool has_vdw_opsl_values(std::string const& residue, std::string const& atom, std::string const& element);