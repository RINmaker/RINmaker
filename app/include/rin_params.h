#pragma once

#include <algorithm>
#include <string>
#include <optional>
#include <variant>
#include <filesystem>

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

    // general stuff
    struct output_file
    { std::filesystem::path value; };

    struct output_directory
    { std::filesystem::path value; };

    enum class illformed_policy_t
    {
        FAIL, SKIP_RES, KEEP_RES, KEEP_ALL
    };

private:
    double _query_dist_hbond = cfg::params::query_dist_hbond;
    double _surface_dist_vdw = cfg::params::surface_dist_vdw;
    double _query_dist_ionic = cfg::params::query_dist_ionic;
    double _query_dist_pipi = cfg::params::query_dist_pipi;
    double _query_dist_pica = cfg::params::query_dist_pica;

    double _query_dist_cmap = cfg::params::query_dist_alpha;

    double _hbond_angle{cfg::params::hbond_angle},
        _pipistack_normal_centre_angle_range{cfg::params::pipistack_normal_centre_angle_range},
        _pipistack_normal_normal_angle_range{cfg::params::pipistack_normal_normal_angle_range},
        _pication_angle{cfg::params::pication_angle};

    int _sequence_separation = cfg::params::seq_sep;

    interaction_type_t _interaction_type = interaction_type_t::NONCOVALENT_BONDS;
    network_policy_t _network_policy = network_policy_t::ALL;
    contact_map_type_t _cmap_type = contact_map_type_t::ALPHA;

    bool _hbond_realistics{true};

    // general stuff
    std::filesystem::path _input{};
    std::variant<output_file, output_directory> _output{};

    bool _skip_water{false}, _no_hydrogen{false};

    illformed_policy_t _illformed{};

    parameters() = default;

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
    auto hbond_angle() const
    { return _hbond_angle; }

    [[nodiscard]]
    auto pipistack_normal_centre_angle_range() const
    { return _pipistack_normal_centre_angle_range; }

    [[nodiscard]]
    auto pipistack_normal_normal_angle_range() const
    { return _pipistack_normal_normal_angle_range; }

    [[nodiscard]]
    auto pication_angle() const
    { return _pication_angle; }

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

    // general stuff
    [[nodiscard]]
    auto const& input() const
    { return _input; }

    [[nodiscard]]
    auto const& output() const
    { return _output; }

    [[nodiscard]]
    auto skip_water() const
    { return _skip_water; }

    [[nodiscard]]
    auto no_hydrogen() const
    { return _no_hydrogen; }

    [[nodiscard]]
    auto illformed_policy() const
    { return _illformed; }
};

struct parameters::configurator final
{
private:
    parameters params;

public:
    [[nodiscard]]
    parameters build() const
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

    configurator& set_hbond_angle(double val)
    {
        params._hbond_angle = val;
        return *this;
    }

    configurator& set_pipistack_normal_centre_angle_range(double val)
    {
        params._pipistack_normal_centre_angle_range = val;
        return *this;
    }

    configurator& set_pipistack_normal_normal_angle_range(double val)
    {
        params._pipistack_normal_normal_angle_range = val;
        return *this;
    }

    configurator& set_pication_angle(double val)
    {
        params._pication_angle = val;
        return *this;
    }

    // general stuff
    configurator& set_input(std::filesystem::path const& path)
    {
        params._input = path;
        return *this;
    }

    configurator& set_output(std::filesystem::path const& path, bool is_directory)
    {
        if (is_directory)
            params._output = rin::parameters::output_directory{path};
        else
            params._output = rin::parameters::output_file{path};

        return *this;
    }

    configurator& set_skip_water(bool skip_water)
    {
        params._skip_water = skip_water;
        return *this;
    }

    configurator& set_no_hydrogen(bool no_hydrogen)
    {
        params._no_hydrogen = no_hydrogen;
        return *this;
    }

    configurator& set_illformed_policy(illformed_policy_t on_malformed)
    {
        params._illformed = on_malformed;
        return *this;
    }
};
}