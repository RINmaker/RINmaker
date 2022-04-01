#pragma once

#include <algorithm>
#include "config.h"

namespace rin
{
struct parameters final
{
public:
    enum class interaction_type_t
    {
        NONCOVALENT_BONDS, ALPHA_BACKBONE, BETA_BACKBONE
    };

    enum class network_policy_t
    {
        ALL, BEST_PER_TYPE, BEST_ONE
    };

private:
    double _query_dist_hbond, _surface_dist_vdw, _query_dist_ionic, _query_dist_pipi, _query_dist_pica;
    double _query_dist_alpha, _query_dist_beta;

    int _sequence_separation;

    interaction_type_t _interaction_type;
    network_policy_t _network_policy;

    bool _hbond_realistics;

    parameters() :
            _query_dist_hbond{cfg::params::hbond_strict},
            _surface_dist_vdw{cfg::params::vdw_strict}, // TODO
            _query_dist_ionic{cfg::params::ionic_strict},
            _query_dist_pipi{cfg::params::pipi_strict},
            _query_dist_pica{cfg::params::pication_strict},
            _query_dist_alpha{cfg::params::ca_distance},
            _query_dist_beta{cfg::params::cb_distance},
            _hbond_realistics{true},
            _sequence_separation{cfg::params::seq_sep},
            _interaction_type{interaction_type_t::NONCOVALENT_BONDS},
            _network_policy{network_policy_t::ALL}
    {}

public:
    struct configurator;

    [[nodiscard]]
    double query_dist_hbond() const
    { return _query_dist_hbond; }

    [[nodiscard]]
    double query_dist_vdw() const
    { return _surface_dist_vdw + 2 * cfg::params::max_vdw_radius; }

    [[nodiscard]]
    double query_dist_ionic() const
    { return _query_dist_ionic; }

    [[nodiscard]]
    double query_dist_pipi() const
    { return _query_dist_hbond; }

    [[nodiscard]]
    double query_dist_pica() const
    { return _query_dist_hbond; }

    [[nodiscard]]
    double query_dist_alpha() const
    { return _query_dist_alpha; }

    [[nodiscard]]
    double query_dist_beta() const
    { return _query_dist_beta; }

    [[nodiscard]]
    int sequence_separation() const
    { return _sequence_separation; }

    [[nodiscard]]
    bool hbond_realistic() const
    { return _hbond_realistics; }

    [[nodiscard]]
    interaction_type_t interaction_type() const
    { return _interaction_type; }

    [[nodiscard]]
    network_policy_t network_policy() const
    { return _network_policy; }
};

struct parameters::configurator final
{
private:
    parameters params;

public:
    parameters build()
    { return params; }

    configurator& set_query_dist_hbond(double val)
    {
        params._query_dist_hbond = std::clamp(val, 0.0, cfg::params::max_limit);
        return *this;
    }

    // TODO
    configurator& set_surface_dist_vdw(double val)
    {
        params._surface_dist_vdw = std::clamp(
                val + 2 * cfg::params::max_vdw_radius, 0.0, cfg::params::max_limit);
        return *this;
    }

    configurator& set_query_dist_ionic(double val)
    {
        params._query_dist_ionic = std::clamp(val, 0.0, cfg::params::max_limit);
        return *this;
    }

    configurator& set_query_dist_pipi(double val)
    {
        params._query_dist_pipi = std::clamp(val, 0.0, cfg::params::max_limit);
        return *this;
    }

    configurator& set_query_dist_pica(double val)
    {
        params._query_dist_pica = std::clamp(val, 0.0, cfg::params::max_limit);
        return *this;
    }

    configurator& set_query_dist_alpha(double val)
    {
        params._query_dist_alpha = std::clamp(val, 0.0, cfg::params::max_limit);
        return *this;
    }

    configurator& set_query_dist_beta(double val)
    {
        params._query_dist_beta = std::clamp(val, 0.0, cfg::params::max_limit);
        return *this;
    }

    configurator& set_interaction_type(interaction_type_t interaction_type)
    {
        params._interaction_type = interaction_type;
        return *this;
    }

    configurator& set_network_policy(network_policy_t network_policy)
    {
        params._network_policy = network_policy;
        return *this;
    }
};
}