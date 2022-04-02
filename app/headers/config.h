#pragma once

#include <string>

namespace constant
{
extern const double ion_ion_k;
extern const double avogadro;
extern const double bohr_radius;
}

namespace cfg
{
namespace log
{
extern const char* main_logger_id;
extern const char* console_logger_id;
extern const char* file_logger_id;
extern const char* default_dirname;
extern const std::string file_logger_filename;
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
extern const double pipistack_normal_normal_angle_range;
extern const double pipistack_normal_centre_angle_range;
extern const double hbond_angle;
extern const double pication_angle;
}

namespace graphml
{
extern char const* null;                // "-999.9"
extern char const* none;                // "None"
extern char const* default_dirname;     // "./outputs"
extern char const* output_file_suffix;  // ".graphml"
}

namespace ver
{
extern const char* app_name;
extern const unsigned int major;
extern const unsigned int minor;
extern const unsigned int rev;
}
}