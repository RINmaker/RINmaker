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
    app.add_option("pdb", pdb_path, "Path to PDB file (.pdb)")
       ->required()
       ->check(CLI::ExistingFile);

    filesystem::path log_dir = cfg::log::default_dirname;
    app.add_option("-l,--log-dir", log_dir, "Log directory")
       ->default_str(cfg::log::default_dirname);

    filesystem::path out_dir = cfg::graphml::default_dirname;
    app.add_option("-o,--out-dir", out_dir, "Output directory")
        ->default_str(cfg::graphml::default_dirname);

    int sequence_separation;
    app.add_option("--seq-sep", sequence_separation, "Minimum sequence separation")
       ->default_val(cfg::params::seq_sep);

    string network_policy;
    app.add_option("--network-policy", network_policy, "Available options: all, multiple, one")
       ->default_val("all")
       ->check(
               [](string const& str)->string
               {
                   return str != "all" && str != "multiple" && str != "one"
                          ? string(R"(must be "all", "multiple" or "one" but you entered ")" + str + "\"")
                          : "";
               });

    string interaction_type;
    app.add_option("--interaction-type", interaction_type, "Available options: closest, ca, cb")
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
    app.add_option("--hydrogen-bond", h_distance, "Query distance for hydrogen bonds")
       ->default_val(cfg::params::query_dist_hbond)
       ->check(positive_check);

    double vdw_distance;
    app.add_option("--vdw-bond", vdw_distance, "Surface distance for vdw bonds")
       ->default_val(cfg::params::surface_dist_vdw)
       ->check(positive_check);

    double ionic_distance;
    app.add_option("--ionic-bond", ionic_distance, "Query distance for ionic bonds")
       ->default_val(cfg::params::query_dist_ionic)
       ->check(positive_check);

    double generic_distance;
    app.add_option("--generic-bond", generic_distance, "Query distance for generic bonds")
       ->default_val(cfg::params::query_dist_alpha)
       ->check(positive_check);

    double pication_distance;
    app.add_option("--pication-bond", pication_distance, "Query distance for cation-pi bonds")
       ->default_val(cfg::params::query_dist_pica)
       ->check(positive_check);

    double pipistack_distance;
    app.add_option("--pipistack-bond", pipistack_distance, "Query distance for pi-pi stackings")
       ->default_val(cfg::params::query_dist_pipi)
       ->check(positive_check);

    bool hbond_realistic_flag = false;
    app.add_flag("--h-bond-realistic", hbond_realistic_flag, "filter hydrogen bonds limiting bond per atom");

    // advanced params
    double hbond_angle;
    app.add_option("--h-bond-angle", hbond_angle, "Angle for hydrogen bonds")
       ->default_val(cfg::params::hbond_angle)
       ->check(positive_check);

    double pication_angle;
    app.add_option("--pication-angle", pication_angle, "Angle for cation-pi bonds")
       ->default_val(cfg::params::pication_angle)
       ->check(positive_check);

    double pipistack_normal_normal_angle_range;
    app.add_option("--pipistack-normal-normal", pipistack_normal_normal_angle_range, "Angle range from normal to normal for pi-pi stackings")
       ->default_val(cfg::params::pipistack_normal_normal_angle_range)
       ->check(positive_check);

    double pipistack_normal_centre_angle_range;
    app.add_option("--pipistack-normal-centre", pipistack_normal_centre_angle_range, "Angle range from normal to centre for pi-pi stackings")
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
        pcfg.set_interaction_type(rin::parameters::interaction_type_t::GENERIC_ALPHA);
    else if (interaction_type == "cb")
        pcfg.set_interaction_type(rin::parameters::interaction_type_t::GENERIC_BETA);
    else if (interaction_type == "closest")
        pcfg.set_interaction_type(rin::parameters::interaction_type_t::NONCOVALENT_BONDS);
    else
        throw runtime_error("incorrect interaction type argument: \"" + interaction_type + "\"");

    fs::create_directory(log_dir);
    log_manager::initialize(log_dir);

    auto params = pcfg.build();

    fs::create_directory(out_dir);
    auto out_path = out_dir / pdb_path.filename().replace_extension(".graphml");

    result = arguments{params, pdb_path, out_path, log_dir};
    return true;
}

std::vector<std::pair<uint32_t, std::string>> read_lines(std::filesystem::path const& file_path)
{
    std::ifstream file;
    file.open(file_path);

    // might throw
    if (!file.is_open())
        throw std::runtime_error("could not open " + file_path.string() + "\n");

    string line;
    uint32_t line_number = 0;
    std::vector<std::pair<uint32_t, std::string>> result;
    while (getline(file, line))
        result.emplace_back(line_number++, line);

    file.close();
    return result;
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
