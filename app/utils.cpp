#include "utils.h"

using namespace std;

string app_full_name() {
    using namespace cfg::ver;
    char buf[500];
    const char* dbg =
#   if NDEBUG
            ""
#   else
            "[DEBUG]"
#   endif
    ;
    const char* os =
#   if _WIN32 && !_WIN64
            "Windows 32 bit"
#   elif _WIN32 && _WIN64
            "Windows 64 bit"
#   elif __APPLE__ || __MACH__ || _MAC
    "OSX"
#   elif __linux__
    "Linux"
#   elif __unix__
"Unix"
#   else
"Unknown OS"
#   endif
    ;

//#   if _WIN32
//    time_t now = time(0);
//    tm local_datetime;
//    localtime_s(&local_datetime , &now);
//    int two_digit_year = (1900 + local_datetime.tm_year) % 100;
//#   elif __linux__
//    time_t now = time(0);
//    tm local_datetime;
//    localtime_r(&now, &local_datetime);
//    int two_digit_year = (1900 + local_datetime.tm_year) % 100;
//#   else
//    int two_digit_year = 22;
//# endif

    string date = __DATE__;
    const char two_digit_year[] = {date[date.size() - 2], date[date.size() - 1], '\0'};

    snprintf(
            buf, sizeof(buf), "%s v%d.%d.%d make_instance %s %s (%s) %s\n(C) 2020-%s Ca' Foscari University of Venice\n", app_name, major, minor, rev, __DATE__, __TIME__, os, dbg
            , two_digit_year); // sprintf for easier format string under C++17 (C++20 would have std::format)
    return buf;
}

bool readArgs(int argc, const char* argv[]) {
    if (argc <= 1) {
        cout << "Use -h or --help for help." << endl;
        return false;
    }

    // CLI stuff
    //

    CLI::App app(app_full_name());
    app.get_formatter()->column_width(48);

    filesystem::path pdb_file;
    app.add_option("path-to-pdb", pdb_file, "PDB (.pdb) file input")
       ->required()
       ->check(CLI::ExistingFile);

    filesystem::path log_dir = cfg::log::default_dirname;
    app.add_option("-l,--log", log_dir, "log directory")
       ->default_str(cfg::log::default_dirname);

    filesystem::path out_file;
    app.add_option("-o,--output", out_file, "graphml (.xml) file output");

    unsigned int seq_sep;
    app.add_option("--seq-sep", seq_sep, "sequence separation")
       ->default_val(cfg::params::seq_sep);

    string bond_control;
    app.add_option("--bond-control", bond_control, "strict or weak")
       ->default_val("strict")
       ->check(
               [](string const& str) {
                   return
                           str != "strict" && str != "weak"
                           ? string(R"(must be "strict" or "weak" but you entered ")" + str + "\"")
                           : "";
               });

    std::string interaction_type;
    app.add_option("--interaction-type", interaction_type, "all, multiple, one")
       ->default_val("all")
       ->check(
               [](string const& str) {
                   if (str != "all" && str != "multiple" && str != "one") {
                       return string("must be \"all\", \"multiple\", \"one\" but you entered \"" + str + "\"");
                   } else {
                       return string();
                   }
               });

    string net_policy;
    app.add_option("--net-policy", net_policy, "closest, ca or cb")
       ->default_val("closest")
       ->check(
               [](string const& str) {
                   if (str != "closest" && str != "ca" && str != "cb") {
                       return string("must be \"closest\", \"ca\" or \"cb\" but you entered \"" + str + "\"");
                   } else {
                       return string();
                   }
               });

    auto positive_check = [](std::string const& str) {
        double val = stod(str);
        if (val <= 0) {
            return std::string("cannot be <= 0");
        } else {
            return std::string();
        }
    };

    double h_distance;
    app.add_option("--h-bond", h_distance, "maximum distance for h bonds")
       ->default_val(cfg::params::hbond_strict)
       ->check(positive_check);

    double vdw_distance;
    app.add_option("--vdw-bond", vdw_distance, "maximum distance for vdw bonds")
       ->default_val(cfg::params::vdw_strict)
       ->check(positive_check);

    double ionic_distance;
    app.add_option("--ionic-bond", ionic_distance, "maximum distance for ionic bonds")
       ->default_val(cfg::params::ionic_strict)
       ->check(positive_check);

    double generic_distance;
    app.add_option("--generic-bond", generic_distance, "maximum distance for generic bonds")
       ->default_val(cfg::params::generic_strict)
       ->check(positive_check);

    double pication_distance;
    app.add_option("--pication-bond", pication_distance, "maximum distance for pication bonds")
       ->default_val(cfg::params::pication_strict)
       ->check(positive_check);

    double pipistack_distance;
    app.add_option("--pipistack-bond", pipistack_distance, "maximum distance for pipistack bonds")
       ->default_val(cfg::params::pipi_strict)
       ->check(positive_check);

    bool force_flag = false;
    app.add_flag("--force", force_flag, "force disable safety guards for distance arguments");

    bool hbond_realistic_flag = false;
    app.add_flag("--h-bond-realistic", hbond_realistic_flag, "filter hydrogen bonds limiting bond per atom");

    // advanced params
    double hbond_angle;
    app.add_option("--h-bond-angle", hbond_angle, "angle for h bonds")
       ->default_val(cfg::params::hbond_angle)
       ->check(positive_check);

    double pication_angle;
    app.add_option("--pication-angle", pication_angle, "angle for pication bonds")
       ->default_val(cfg::params::pication_angle)
       ->check(positive_check);

    double pipistack_normal_normal_angle_range;
    app.add_option("--pipistack-normal-normal", pipistack_normal_normal_angle_range, "angle range from normal to normal for pipistack bonds")
       ->default_val(cfg::params::pipistack_normal_normal_angle_range)
       ->check(positive_check);

    double pipistack_normal_centre_angle_range;
    app.add_option("--pipistack-normal-centre", pipistack_normal_centre_angle_range, "angle range from normal to centre for pipistack bonds")
       ->default_val(cfg::params::pipistack_normal_centre_angle_range)
       ->check(positive_check);

    // CLI parser
    CLI11_PARSE(app, argc, argv);

    // init log
    //
    parameters::set_log_path(log_dir);

    parameters::set_pdb_path(pdb_file);
    parameters::set_out_path(out_file);

    parameters::set_net_policy(net_policy);
    parameters::set_bond_control(bond_control);
    parameters::set_interaction_type(interaction_type);

    parameters::set_forcing(force_flag);
    parameters::set_hbond_realistic(hbond_realistic_flag);

    parameters::set_hbond_distance(h_distance);
    parameters::set_vdw_distance(vdw_distance);
    parameters::set_ionic_distance(ionic_distance);
    parameters::set_generic_distance(generic_distance);
    parameters::set_pication_distance(pication_distance);
    parameters::set_pipistack_distance(pipistack_distance);

    parameters::set_seq_sep(seq_sep);

    std::filesystem::create_directory(log_dir);
    log_manager::initialize(log_dir);

    return true;
}
