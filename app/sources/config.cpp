#include <string>
#include "config.h"
#include <limits>

using std::string;

namespace constant
{
const double ion_ion_k = 33.4450992;
const double avogadro = 6.02214076E23;
const double bohr_radius = 5.2917721E-11;//m
const double pipi_a = -0.5274;
const double pipi_b = 25.6290;
const double pipi_c = -25.653;
}

namespace cfg
{
namespace ver
{
const char* app_name = "RING";
const unsigned int major = 0;
const unsigned int minor = 1;
const unsigned int rev = 3;
}

namespace log
{
char const* const main_logger_id = "main";
char const* const console_logger_id = "console";
char const* const file_logger_id = "file";
char const* const default_file_logger_filename = "./main.txt";
}

namespace params
{
const double query_dist_hbond = 3.5;
const double surface_dist_vdw = 0.5;
const double query_dist_ionic = 4.0;
const double query_dist_pipi = 6.5;
const double query_dist_pica = 5.0;

const double query_dist_alpha = 6.0;
const double query_dist_beta = 6.0;

const int seq_sep = 3;

const double max_limit = 20.0;

const double max_vdw_radius = 1.90;
const double max_pipi_atom_atom_distance = 4.5;

const double pipistack_normal_normal_angle_range = 30;
const double pipistack_normal_centre_angle_range = 60;

const double hbond_angle = 63.;
const double pication_angle = 45.;
}

namespace graphml
{
char const* null = "-999.9";
char const* none = "None";
char const* default_dirname = ".";
char const* default_filename = "./network.graphml";
char const* output_file_suffix = ".graphml";
}
}