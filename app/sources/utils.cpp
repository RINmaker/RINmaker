#include "utils.h"

using namespace std;

namespace fs = std::filesystem;

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

    string date = __DATE__; //Date at compile time
    const char two_digit_year[] = {date[date.size() - 2], date[date.size() - 1], '\0'};
    snprintf(
            buf, sizeof(buf), "%s v%d.%d.%d build %s %s (%s) %s\n(C) 2020-%s Ca' Foscari University of Venice\n", app_name, major, minor, rev, __DATE__, __TIME__, os, dbg
            , two_digit_year); // sprintf for easier format string under C++17 (C++20 would have std::format)
    return buf;
}

optional<arguments> read_args(int argc, const char* argv[])
{
    if (argc <= 1)
    {
        cout << "Use -h or --help for help." << endl;
        return nullopt;
    }

    CLI::App app(app_full_name());
    app.set_help_all_flag("-H,--help-expanded", "Print this help message (expanded) and exit");
    app.get_formatter()->column_width(64);

    filesystem::path pdb_path;
    app.add_option("-i, --input", pdb_path, "Path to PDB file (.pdb)")
        ->required()
        ->check(CLI::ExistingFile);

    filesystem::path log_dir = cfg::log::default_dirname;
    app.add_option("-l,--log-directory", log_dir, "Log directory")
        ->default_str(cfg::log::default_dirname);

    filesystem::path out_dir = cfg::graphml::default_dirname;
    app.add_option("-o,--out-directory", out_dir, "Output directory")
        ->default_str(cfg::graphml::default_dirname);

    uint32_t sequence_separation;
    app.add_option("-s,--sequence-separation", sequence_separation, "Minimum sequence separation")
        ->default_val(cfg::params::seq_sep)
        ->check(CLI::PositiveNumber);

    // rin subcommand
    auto rin_app = app.add_subcommand(
            "rin", "Compute the residue interaction network");

    auto network_policy = rin::parameters::network_policy_t::ALL;
    std::map<std::string, rin::parameters::network_policy_t> netp_map{
            {"all", rin::parameters::network_policy_t::ALL},
            {"one", rin::parameters::network_policy_t::BEST_ONE},
            {"multiple", rin::parameters::network_policy_t::BEST_PER_TYPE}};

    rin_app->add_option(
            "--policy", network_policy, "Affects which edges are kept per pair of aminoacids")
            ->transform(
                    CLI::CheckedTransformer(netp_map, CLI::ignore_case).description(
                            CLI::detail::generate_map(CLI::detail::smart_deref(netp_map), true)))
            ->default_val("all");

    double h_distance;
    rin_app->add_option("--hydrogen-bond", h_distance, "Query distance for hydrogen bonds")
       ->default_val(cfg::params::query_dist_hbond)
       ->check(CLI::PositiveNumber);

    double vdw_distance;
    rin_app->add_option("--vdw-bond", vdw_distance, "Surface distance for vdw bonds")
       ->default_val(cfg::params::surface_dist_vdw)
       ->check(CLI::PositiveNumber);

    double ionic_distance;
    rin_app->add_option("--ionic-bond", ionic_distance, "Query distance for ionic bonds")
       ->default_val(cfg::params::query_dist_ionic)
       ->check(CLI::PositiveNumber);

    double pication_distance;
    rin_app->add_option("--pication-bond", pication_distance, "Query distance for cation-pi bonds")
       ->default_val(cfg::params::query_dist_pica)
       ->check(CLI::PositiveNumber);

    double pipistack_distance;
    rin_app->add_option("--pipistack-bond", pipistack_distance, "Query distance for pi-pi stackings")
       ->default_val(cfg::params::query_dist_pipi)
       ->check(CLI::PositiveNumber);

    bool hbond_realistic_flag = false;
    rin_app->add_flag("--h-bond-realistic", hbond_realistic_flag, "filter hydrogen bonds limiting bond per atom");

    // advanced params
    double hbond_angle;
    rin_app->add_option("--h-bond-angle", hbond_angle, "Angle for hydrogen bonds")
       ->default_val(cfg::params::hbond_angle)
       ->check(CLI::PositiveNumber);

    double pication_angle;
    rin_app->add_option("--pication-angle", pication_angle, "Angle for cation-pi bonds")
        ->default_val(cfg::params::pication_angle)
        ->check(CLI::PositiveNumber);

    double pipistack_normal_normal_angle_range;
    rin_app->add_option("--pipistack-normal-normal", pipistack_normal_normal_angle_range, "Angle range from normal to normal for pi-pi stackings")
       ->default_val(cfg::params::pipistack_normal_normal_angle_range)
       ->check(CLI::PositiveNumber);

    double pipistack_normal_centre_angle_range;
    rin_app->add_option("--pipistack-normal-centre", pipistack_normal_centre_angle_range, "Angle range from normal to centre for pi-pi stackings")
       ->default_val(cfg::params::pipistack_normal_centre_angle_range)
       ->check(CLI::PositiveNumber);

    // contact map subcommand
    auto cmap_app = app.add_subcommand(
            "cmap", "Compute the contact map of the protein");

    auto cmap_type = rin::parameters::contact_map_type_t::ALPHA;
    std::map<std::string, rin::parameters::contact_map_type_t> cmt_map{
            {"ca", rin::parameters::contact_map_type_t::ALPHA},
            {"cb", rin::parameters::contact_map_type_t::BETA}};

    cmap_app->add_option(
                    "--type", cmap_type, "Type of contact map (alpha/beta carbon)")
            ->transform(
                    CLI::CheckedTransformer(cmt_map, CLI::ignore_case).description(
                            CLI::detail::generate_map(CLI::detail::smart_deref(cmt_map), true)))
            ->default_val("ca");

    double generic_distance;
    cmap_app->add_option("--distance", generic_distance, "Query distance between alpha/beta carbons")
            ->default_val(cfg::params::query_dist_alpha)
            ->check(CLI::PositiveNumber);

    // rin XOR cmap REQUIRED
    app.require_subcommand(1, 1);

    //CLI11_PARSE(app, argc, argv);
    try {
        app.parse(argc, argv);
    } catch(const CLI::ParseError &e) {
        app.exit(e);
        return nullopt;
    }

    auto pcfg = rin::parameters::configurator()
            .set_query_dist_cmap(generic_distance)

            .set_query_dist_hbond(h_distance)
            .set_surface_dist_vdw(vdw_distance)
            .set_query_dist_ionic(ionic_distance)
            .set_query_dist_pica(pication_distance)
            .set_query_dist_pipi(pipistack_distance)

            .set_sequence_separation(sequence_separation)
            .set_hbond_realistic(hbond_realistic_flag)

            .set_network_policy(network_policy)
            .set_cmap_type(cmap_type);

    if(rin_app->parsed())
        pcfg.set_interaction_type(rin::parameters::interaction_type_t::NONCOVALENT_BONDS);

    else //if(cmap_app->parsed())
        pcfg.set_interaction_type(rin::parameters::interaction_type_t::CONTACT_MAP);

    fs::create_directory(log_dir);
    log_manager::initialize(log_dir);

    auto params = pcfg.build();

    fs::create_directory(out_dir);

    return arguments{params, pdb_path, out_dir, log_dir};
}

string joinStrings(std::vector<std::string> const& values, string const& delimiter)
{
    string out;
    for (string const& value: values)
    {
        if (!out.empty())
            out += delimiter;
        out += value;
    }

    return out;
}
