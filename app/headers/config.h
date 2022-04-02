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
extern const char* output_file_suffix;

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
extern const double hbond_strict;
extern const double vdw_strict;
extern const double ionic_strict;
extern const double pipi_strict;
extern const double pication_strict;

extern const double ca_distance;
extern const double cb_distance;

extern const double generic_strict;

extern char const* _bond_control;
extern char const* net_policy;
extern char const* interaction_type;

extern const int seq_sep;

extern const double weak_powering;
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
extern char const* null;            // "-999.9"
extern char const* none;            // "None"
extern char const* default_dirname; // "./outputs"
}

namespace ver
{
extern const char* app_name;
extern const unsigned int major;
extern const unsigned int minor;
extern const unsigned int rev;
}
}