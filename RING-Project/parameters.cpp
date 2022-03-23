#include "parameters.h"
#include "config.h"

#include <sstream>

std::string __to_string(parameters::policy p)
{
    switch (p)
    {
    case parameters::policy::CLOSEST:
        return "\"closest\"";

    case parameters::policy::CA:
        return "\"ca\"";

    case parameters::policy::CB:
        return "\"cb\"";

    default:
        return "\"unknown\"";
    }
}

std::string __to_string(parameters::control c)
{
    switch (c)
    {
    case parameters::control::STRICT_BOND:
        return "\"strict\"";

    case parameters::control::WEAK_BOND:
        return "\"weak\"";

    default:
        return "\"unknown\"";
    }
}

std::string __to_string(parameters::interaction_type i)
{
    switch (i)
    {
    case parameters::interaction_type::ALL:
        return "\"all\"";

    case parameters::interaction_type::MULTIPLE:
        return "\"multiple\"";

    case parameters::interaction_type::ONE:
        return "\"one\"";

    default:
        return "\"unknown\"";
    }
}

std::string parameters::pretty()
{
    std::ostringstream strs;
    
    strs << "{";

    policy const pc = get_net_policy();
    strs << "net_policy:" << __to_string(pc) << ", ";
    switch (pc)
    {
    case policy::CLOSEST:
        strs
            << "h_bond:"        << get_distance_h()         << ", "
            << "vdw_bond:"      << get_distance_vdw()       << ", "
            << "ionic_bond:"    << get_distance_ionic()     << ", "
            << "pication:"      << get_distance_pication()  << ", "
            << "pipistack:"     << get_distance_pipistack() << ", ";
        break;

    case policy::CA:
    case policy::CB:
        strs << "generic_bond:" << get_distance_generic()   << ", ";
        break;

    default:
        strs << "\"unknown\",";
        break;
    }

    strs
        << "seq_sep:"           << get_seq_sep()                        << ", "
        << "bond_control:"      << __to_string(get_bond_control())      << ", "
        << "interaction_type:"  << __to_string(get_interaction_type());
    
    strs << "}";

    return strs.str();
}