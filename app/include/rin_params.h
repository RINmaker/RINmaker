#pragma once

#include <algorithm>
#include <string>
#include "config.h"

namespace rin
{
struct parameters final
{
public:
    enum class interaction_type_t
    {
        NONCOVALENT_BONDS, CONTACT_MAP
    };

    enum class network_policy_t
    {
        ALL, BEST_PER_TYPE, BEST_ONE
    };

    enum class contact_map_type_t
    {
        ALPHA, BETA
    };

private:
    double _query_dist_hbond, _surface_dist_vdw, _query_dist_ionic, _query_dist_pipi, _query_dist_pica;
    double _query_dist_cmap;

    int _sequence_separation;

    interaction_type_t _interaction_type;
    network_policy_t _network_policy;
    contact_map_type_t _cmap_type;

    bool _hbond_realistics;

    parameters() :
            _query_dist_hbond{cfg::params::query_dist_hbond},
            _surface_dist_vdw{cfg::params::surface_dist_vdw},
            _query_dist_ionic{cfg::params::query_dist_ionic},
            _query_dist_pipi{cfg::params::query_dist_pipi},
            _query_dist_pica{cfg::params::query_dist_pica},
            _query_dist_cmap{cfg::params::query_dist_alpha},
            _hbond_realistics{true},
            _sequence_separation{cfg::params::seq_sep},
            _interaction_type{interaction_type_t::NONCOVALENT_BONDS},
            _network_policy{network_policy_t::ALL},
            _cmap_type{contact_map_type_t::ALPHA}
    {}

public:
    struct configurator;

    [[nodiscard]]
    double query_dist_hbond() const
    { return _query_dist_hbond; }

    [[nodiscard]]
    double query_dist_vdw() const
    { return std::clamp(_surface_dist_vdw + 2 * cfg::params::max_vdw_radius, 0.0, cfg::params::max_limit); }

    [[nodiscard]]
    double surface_dist_vdw() const
    { return _surface_dist_vdw; }

    [[nodiscard]]
    double query_dist_ionic() const
    { return _query_dist_ionic; }

    [[nodiscard]]
    double query_dist_pipi() const
    { return _query_dist_pipi; }

    [[nodiscard]]
    double query_dist_pica() const
    { return _query_dist_pica; }

    [[nodiscard]]
    double query_dist_cmap() const
    { return _query_dist_cmap; }

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

    [[nodiscard]]
    contact_map_type_t cmap_type() const
    { return _cmap_type; }

    [[nodiscard]]
    std::string pretty() const;
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

    configurator& set_surface_dist_vdw(double val)
    {
        params._surface_dist_vdw = val;
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

    configurator& set_query_dist_cmap(double val)
    {
        params._query_dist_cmap = std::clamp(val, 0.0, cfg::params::max_limit);
        return *this;
    }

    configurator& set_interaction_type(interaction_type_t interaction_type)
    {
        params._interaction_type = interaction_type;
        return *this;
    }

    configurator& set_cmap_type(contact_map_type_t cmap_type)
    {
        params._cmap_type = cmap_type;
        return *this;
    }

    configurator& set_network_policy(network_policy_t network_policy)
    {
        params._network_policy = network_policy;
        return *this;
    }

    configurator& set_sequence_separation(int val)
    {
        params._sequence_separation = std::max(cfg::params::seq_sep, val);
        return *this;
    }

    configurator& set_hbond_realistic(bool val)
    {
        params._hbond_realistics = val;
        return *this;
    }
};
}