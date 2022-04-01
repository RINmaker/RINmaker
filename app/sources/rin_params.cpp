#include "rin_params.h"
#include "config.h"

#include <sstream>

using netp = rin::parameters::network_policy_t;
using intt = rin::parameters::interaction_type_t;
using std::string;

string to_string(intt i)
{
    switch (i)
    {
    case intt::NONCOVALENT_BONDS:
        return "\"closest\"";

    case intt::ALPHA_BACKBONE:
        return "\"ca\"";

    case intt::BETA_BACKBONE:
        return "\"cb\"";

    default:
        return "\"unknown\"";
    }
}

string to_string(netp p)
{
    switch (p)
    {
    case netp::ALL:
        return "\"all\"";

    case netp::BEST_PER_TYPE:
        return "\"multiple\"";

    case netp::BEST_ONE:
        return "\"one\"";

    default:
        return "\"unknown\"";
    }
}

string rin::parameters::pretty() const
{
    std::ostringstream strs;

    strs << "{";

    auto pc = interaction_type();
    strs << "net_policy:" << to_string(pc) << ", ";
    switch (pc)
    {
    case intt::NONCOVALENT_BONDS:
        strs
                << "h_bond:" << query_dist_hbond() << ", "
                << "vdw_bond:" << query_dist_vdw() << ", "
                << "ionic_bond:" << query_dist_ionic() << ", "
                << "pication:" << query_dist_pica() << ", "
                << "pipistack:" << query_dist_pipi() << ", ";
        break;

    case intt::ALPHA_BACKBONE:
        strs << "generic_bond:" << query_dist_alpha() << ", ";
        break;

    case intt::BETA_BACKBONE:
        strs << "generic_bond:" << query_dist_beta() << ", ";
        break;
    }

    strs << "seq_sep:" << sequence_separation() << ", " << "interaction_type:" << to_string(network_policy()) << "}";

    return strs.str();
}