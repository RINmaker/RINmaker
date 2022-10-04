#pragma once

#include <string>

namespace constant
{
extern const double ion_ion_k;
extern const double avogadro;
extern const double bohr_radius;
extern const double pipi_a;
extern const double pipi_b;
extern const double pipi_c;
}

namespace cfg
{
extern char const* const monomer_lib_dir;

namespace log
{
extern char const* const main_logger_id;
extern char const* const console_logger_id;
extern char const* const file_logger_id;
extern char const* const default_file_logger_filename;
}

namespace params
{
extern const double query_dist_hbond;
extern const double surface_dist_vdw;
extern const double query_dist_ionic;
extern const double query_dist_pipi;
extern const double query_dist_pica;

extern const double query_dist_alpha;
extern const double query_dist_beta;

extern const int seq_sep;

extern const double max_limit;

extern const double max_vdw_radius;
extern const double max_pipi_atom_atom_distance;

// advanced parameters for deep testing
extern double const pipistack_normal_normal_angle_range;
extern double const pipistack_normal_centre_angle_range;
extern double const hbond_angle;
extern double const pication_angle;
}

namespace graphml
{
extern char const* const null;                // "-999.9"
extern char const* const none;                // "None"
extern char const* const default_dirname;     // "."
extern char const* const default_filename;    // "./network.graphml
extern char const* const output_file_suffix;  // ".graphml"
}

namespace ver
{
extern char const* const app_name;
extern unsigned int const major;
extern unsigned int const minor;
extern unsigned int const rev;
}
}