#include "rin_params.h"
#include "config.h"

#include <sstream>

using std::string;

string to_string(rin::parameters::network_policy_t network_policy)
{
    switch (network_policy)
    {
    case rin::parameters::network_policy_t::ALL:
        return "all";

    case rin::parameters::network_policy_t::BEST_PER_TYPE:
        return "multiple";

    case rin::parameters::network_policy_t::BEST_ONE:
        return "one";

    default:
        // unreachable
        return "err";
    }
}


string to_string(rin::parameters::interaction_type_t interaction_type)
{
    switch (interaction_type)
    {
        case rin::parameters::interaction_type_t::NONCOVALENT_BONDS:
            return "closest";

        case rin::parameters::interaction_type_t::CONTACT_MAP:
            return "cmap";

        default:
            // unreachable
            return "err";
    }
}

string to_string(rin::parameters::contact_map_type_t cmap_type)
{
    switch (cmap_type)
    {
    case rin::parameters::contact_map_type_t::ALPHA:
        return "alpha";

    case rin::parameters::contact_map_type_t::BETA:
        return "beta";

    default:
        // unreachable
        return "err";
    }
}

string rin::parameters::pretty() const
{
    std::ostringstream strs;

    strs << "{";

    auto in = interaction_type();
    strs << "interaction_type:" << to_string(in) << ", ";
    switch (in)
    {
    case rin::parameters::interaction_type_t::NONCOVALENT_BONDS:
        strs
                << "network_policy:" << to_string(network_policy()) << ","
                << "h_bond:" << query_dist_hbond() << ", "
                << "vdw_bond:" << query_dist_vdw() << ", "
                << "ionic_bond:" << query_dist_ionic() << ", "
                << "pication:" << query_dist_pica() << ", "
                << "pipistack:" << query_dist_pipi() << ", ";
        break;

    case rin::parameters::interaction_type_t::CONTACT_MAP:
        strs
            << "contact_map_type:" << to_string(cmap_type()) << ", "
            << "generic_bond:" << query_dist_cmap() << ", ";
        break;
    }

    strs << "seq_sep:" << sequence_separation() << "}";

    return strs.str();
}