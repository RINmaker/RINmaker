#include "rin_params.h"
#include "config.h"

#include <sstream>

using netp = rin::parameters::network_policy_t;
using intt = rin::parameters::interaction_type_t;
using std::string;

string rin::parameters::to_string(netp np)
{
    switch (np)
    {
    case netp::ALL:
        return "all";

    case netp::BEST_PER_TYPE:
        return "multiple";

    case netp::BEST_ONE:
        return "one";

    default:
        return "unknown";
    }
}

string rin::parameters::to_string(intt it)
{
    switch (it)
    {
    case intt::NONCOVALENT_BONDS:
        return "non_covalent";

    case intt::GENERIC_ALPHA:
        return "generic_alpha";

    case intt::GENERIC_BETA:
        return "generic_beta";

    default:
        return "unknown";
    }
}

string rin::parameters::pretty() const
{
    std::ostringstream strs;

    strs << "{";

    auto pc = interaction_type();
    strs << "interaction_type:" << rin::parameters::to_string(pc) << ", ";
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

    case intt::GENERIC_ALPHA:
        strs << "generic_bond:" << query_dist_alpha() << ", ";
        break;

    case intt::GENERIC_BETA:
        strs << "generic_bond:" << query_dist_beta() << ", ";
        break;
    }

    strs << "seq_sep:" << sequence_separation() << ", " << "network_policy:" << rin::parameters::to_string(network_policy()) << "}";

    return strs.str();
}