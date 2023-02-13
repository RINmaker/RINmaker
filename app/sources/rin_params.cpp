#include "rin_params.h"
#include "config.h"

#include <sstream>

using std::string;

string to_string(rin::parameters::network_policy_t network_policy)
{
    string ret{};
    switch (network_policy)
    {
    case rin::parameters::network_policy_t::ALL:
        ret = "\"all\"";
        break;
    case rin::parameters::network_policy_t::BEST_PER_TYPE:
        ret = "\"multiple\"";
        break;
    case rin::parameters::network_policy_t::BEST_ONE:
        ret = "\"one\"";
        break;
    }
    return ret;
}

string to_string(rin::parameters::contact_map_type_t cmap_type)
{
    string ret{};
    switch (cmap_type)
    {
    case rin::parameters::contact_map_type_t::ALPHA:
        ret = "\"ca\"";
        break;
    case rin::parameters::contact_map_type_t::BETA:
        ret = "\"cb\"";
        break;
    }
    return ret;
}

string to_string(rin::parameters::illformed_policy_t illformed_policy)
{
    string ret{};
    switch (illformed_policy)
    {
    case rin::parameters::illformed_policy_t::FAIL:
        ret = "\"fail\"";
        break;
    case rin::parameters::illformed_policy_t::SKIP_RES:
        ret = "\"sres\"";
        break;
    case rin::parameters::illformed_policy_t::KEEP_RES:
        ret = "\"kres\"";
        break;
    case rin::parameters::illformed_policy_t::KEEP_ALL:
        ret = "\"kall\"";
        break;
    }
    return ret;
}

string rin::parameters::serialize_rin() const
{
    return (std::ostringstream{}
        << "{"
        << "\"policy\": " << to_string(network_policy()) << ", "
        << "\"hydrogen-bond\": " << query_dist_hbond() << ", "
        << "\"vdw-bond\": " << surface_dist_vdw() << ", "
        << "\"ionic-bond\": " << query_dist_ionic() << ", "
        << "\"pication-bond\": " << query_dist_pica() << ", "
        << "\"pipistack-bond\": " << query_dist_pipi() << ", "
        << "\"h-bond-realistic\": " << (hbond_realistic() ? "true" : "false") << ", "
        << "\"h-bond-angle\": " << hbond_angle() << ", "
        << "\"pication-angle\": " << pication_angle() << ", "
        << "\"pipistack-normal-normal\": " << pipistack_normal_normal_angle_range() << ", "
        << "\"pipistack-normal-centre\": " << pipistack_normal_centre_angle_range() << "}"
    ).str();
}

string rin::parameters::serialize_cmap() const
{
    return (std::ostringstream{}
        << "{"
        << "\"type\": " << to_string(cmap_type()) << ","
        << "\"distance\": " << query_dist_cmap() << "}"
    ).str();
}

string rin::parameters::pretty() const
{
    std::ostringstream strs{};
    strs
        << "{"
        << "\"input\":" << '\"' << input().string() << '\"' << ", ";

    if(auto out = output(); std::holds_alternative<rin::parameters::output_directory> (out))
    {
        strs
            << R"("output": ")" << std::get<rin::parameters::output_directory>(out).value.string() << "\", "
            << "\"-d\": true, " ;
    }
    else
    {
        strs
            << R"("output": ")" << std::get<rin::parameters::output_directory>(out).value.string() << "\", "
            << "\"-d\": false, ";
    }

    strs
        << "\"no-hydrogen\": " << (no_hydrogen() ? "true" : "false") << ", "
        << "\"keep-water\": " << (skip_water() ? "false" : "true") << ", "
        << "\"sequence-separation\": " << sequence_separation() << ", "
        << "\"illformed\": " << to_string(illformed_policy()) << ", ";

    switch (interaction_type())
    {
    case rin::parameters::interaction_type_t::NONCOVALENT_BONDS:
        strs << "\"rin\": " << serialize_rin();
        break;

    case rin::parameters::interaction_type_t::CONTACT_MAP:
        strs << "\"cmap\": " << serialize_cmap();
        break;
    }

    strs << "}";
    return strs.str();
}