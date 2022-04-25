#include <string>
#include "config.h"
#include <limits>

using std::string;

namespace constant
{
const double ion_ion_k = 33.4450992;
const double avogadro = 6.02214076E23;
const double bohr_radius = 5.2917721E-11;//m
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
const char* main_logger_id = "main";
const char* console_logger_id = "console";
const char* file_logger_id = "file";

const char* default_dirname = "./logs";
const string file_logger_filename = string(main_logger_id) + ".txt";
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
char const* output_file_suffix = ".graphml";
}
}