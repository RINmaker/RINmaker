#include "utils.h"

using namespace std;

string app_full_name()
{
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
            buf, sizeof(buf), "%s v%d.%d.%d build %s %s (%s) %s\n(C) 2020-%s Ca' Foscari University of Venice\n", app_name, major, minor, rev, __DATE__, __TIME__, os, dbg
            , two_digit_year); // sprintf for easier format string under C++17 (C++20 would have std::format)
    return buf;
}

bool read_args(int argc, const char* argv[], optional<arguments>& result)
{
    result = nullopt;
    if (argc <= 1)
    {
        cout << "Use -h or --help for help." << endl;
        return false;
    }

    CLI::App app(app_full_name());
    app.get_formatter()->column_width(48);

    filesystem::path pdb_path;
    app.add_option("path-to-pdb", pdb_path, "path to PDB (.pdb) input")
       ->required()
       ->check(CLI::ExistingFile);

    filesystem::path log_dir = cfg::log::default_dirname;
    app.add_option("-l,--log", log_dir, "log directory")
       ->default_str(cfg::log::default_dirname);

    filesystem::path out_path;
    app.add_option("-o,--output", out_path, "path to graphml (.xml) output");

    int sequence_separation;
    app.add_option("--seq-sep", sequence_separation, "sequence separation")
       ->default_val(cfg::params::seq_sep);

    string network_policy;
    app.add_option("--network-policy", network_policy, "all, multiple, one")
       ->default_val("all")
       ->check(
               [](string const& str)->string
               {
                   return str != "all" && str != "multiple" && str != "one"
                          ? string(R"(must be "all", "multiple", "one" but you entered ")" + str + "\"")
                          : "";
               });

    string interaction_type;
    app.add_option("--interaction-type", interaction_type, "closest, ca or cb")
       ->default_val("closest")
       ->check(
               [](string const& str)->string
               {
                   return str != "closest" && str != "ca" && str != "cb"
                          ? string(R"(must be "closest", "ca" or "cb" but you entered ")" + str + "\"")
                          : "";
               });

    auto positive_check = [](string const& str)->string
    {
        double val = stod(str);
        return val <= 0 ? "must be > 0" : "";
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

    CLI11_PARSE(app, argc, argv);

    auto pcfg = rin::parameters::configurator()
            .set_query_dist_alpha(generic_distance)
            .set_query_dist_beta(generic_distance)

            .set_query_dist_hbond(h_distance)
            .set_surface_dist_vdw(vdw_distance)
            .set_query_dist_ionic(ionic_distance)
            .set_query_dist_pica(pication_distance)
            .set_query_dist_pipi(pipistack_distance)

            .set_sequence_separation(sequence_separation)
            .set_hbond_realistic(hbond_realistic_flag);

    if (network_policy == "all")
        pcfg.set_network_policy(rin::parameters::network_policy_t::ALL);
    else if (network_policy == "multiple")
        pcfg.set_network_policy(rin::parameters::network_policy_t::BEST_PER_TYPE);
    else if (network_policy == "one")
        pcfg.set_network_policy(rin::parameters::network_policy_t::BEST_ONE);
    else
        throw std::runtime_error("incorret network policy argument: \"" + network_policy + "\"");

    if (interaction_type == "ca")
        pcfg.set_interaction_type(rin::parameters::interaction_type_t::ALPHA_BACKBONE);
    else if (interaction_type == "cb")
        pcfg.set_interaction_type(rin::parameters::interaction_type_t::BETA_BACKBONE);
    else if (interaction_type == "closest")
        pcfg.set_interaction_type(rin::parameters::interaction_type_t::NONCOVALENT_BONDS);
    else
        throw std::runtime_error("incorrect interaction type argument: \"" + interaction_type + "\"");

    std::filesystem::create_directory(log_dir);
    log_manager::initialize(log_dir);

    auto params = pcfg.build();
    rin::parameters::global::instance().set(params);

    result = arguments{params, pdb_path, out_path, log_dir};
    return true;
}
