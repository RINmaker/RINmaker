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
		const double hbond_strict = 3.5;
		const double vdw_strict = 0.5;
		const double ionic_strict = 4.0;
		const double pipi_strict = 6.5;
		const double pication_strict = 5.0;

		const double ca_distance = 6.0;
		const double cb_distance = 6.0;

		const double generic_strict = 6.0;

		char const* _bond_control = "strict";
		char const* net_policy = "closest";
		char const* interaction_type = "all";

		const int seq_sep = 3;

		const double weak_powering = 1.0;
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
		char const* null = "-999.9";    // "-999.9"
		char const* none = "None";      // "None"
	}

	const char* output_file_suffix = ".xml";
}