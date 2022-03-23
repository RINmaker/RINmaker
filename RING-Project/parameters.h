#pragma once

#include <string>
#include <memory>
#include <filesystem>

#include "log_manager.h"
#include "config.h"

// this is a singleton that reunites the runtime configuration of the program.
// its purpose is to check and store the parameters that are input by the user
// and to be read from throughout the rest of the codebase.
//
class parameters final
{
public:
    // this is an enum that dictates the types of interactions
    // that will be sought after at runtime.
    //
    enum class policy
    {
        // ss-bond,h-bond,vdw-bond,ionic-bond,pipi-stacking,cation-pi
        CLOSEST,

        // alpha-carbon
        CA,

        // beta-carbon
        CB,
    };

    // this allows for zero or some tolerance on the upper bound of the query distance.
    //
    enum class control
    {
        // no tolerance, user input is the actual limit
        STRICT_BOND,

        // adds x angstrom of tolerance (x is fixed and defined in config.cpp)
        WEAK_BOND,
    };

    // specifies the types of bonds that will be kept in the 
    // generated rin between each pair of interacting aminoacids.
    //
    enum class interaction_type
    {
        // all valid bonds are kept
        ALL,

        // only the best bond per type is kept
        MULTIPLE,

        // only the best bond is kept
        ONE,
    };

private:
    // query distances
    //
    double _hbond_dist;
    double _vdw_dist;
    double _ionic_dist;
    double _generic_dist;
    double _pication_dist;
    double _pipistack_dist;

    // sequence separation is the minimum allowed positional distance between two
    // aminoacids on the same chain in order to be valid interaction candidates
    //
    int _seq_sep;

    // if !_forcing then the maximum query distance that user can specify
    // is the one defined in config.h/cpp, else there are no limits (w: memory consumption)
    //
    bool _forcing;

    // these paths can be provided by user or left as specified in config
    //
    std::filesystem::path _pdb_path, _output_path, _log_path;

    bool _hbond_realistic;

    policy _net_policy;

    control _bond_control;

    interaction_type _interaction_type;

    [[nodiscard]]
    static parameters& instance()
    {
        static parameters params;
        return params;
    }

    // defaults are specified in config.h/cpp
    //
    parameters()
        : _hbond_dist(cfg::params::hbond_strict)
        , _vdw_dist(cfg::params::vdw_strict)
        , _ionic_dist(cfg::params::ionic_strict)
        , _generic_dist(cfg::params::generic_strict)
        , _pication_dist(cfg::params::pication_strict)
        , _pipistack_dist(cfg::params::pipi_strict)

        , _net_policy(policy::CLOSEST)                  // TODO config
        , _bond_control(control::STRICT_BOND)           // TODO config
        , _interaction_type(interaction_type::MULTIPLE) // TODO config

        , _seq_sep(cfg::params::seq_sep)
        , _forcing(false)                               // TODO config
        , _hbond_realistic(true)

        , _pdb_path("./.pdb")                           // TODO config
        , _output_path("./generated/rin.xml")           // TODO config
        , _log_path("./generated/logs")                 // TODO config
    {}

    // returns the corrected distance value, after applying 
    // the limits defined by config.h/cpp and user input
    //
    [[nodiscard]]
    double __apply_limit(double value) const
    {
        if (_bond_control == control::WEAK_BOND)
        {
            value += cfg::params::weak_powering;
        }

        if (!_forcing && value > cfg::params::max_limit)
        {
            value = cfg::params::max_limit;
        }

        return value;
    }

public:
    // sets the corrected h-bond distance
    //
    static void set_hbond_distance(double d)
    {
        instance()._hbond_dist = instance().__apply_limit(d);
    }

    // sets the corrected vdw-bond distance
    //
    static void set_vdw_distance(double d)
    {
        instance()._vdw_dist = instance().__apply_limit(d);
    }

    // sets the corrected ionic-bond distance
    //
    static void set_ionic_distance(double d)
    {
        instance()._ionic_dist = instance().__apply_limit(d);
    }

    // sets the corrected generic-bond distance
    //
    static void set_generic_distance(double d)
    {
        instance()._generic_dist = instance().__apply_limit(d);
    }

    // sets the corrected cation-pi distance
    //
    static void set_pication_distance(double d)
    {
        instance()._pication_dist = instance().__apply_limit(d);
    }

    // sets the corrected pipi-stacking distance
    //
    static void set_pipistack_distance(double d)
    {
        instance()._pipistack_dist = instance().__apply_limit(d);
    }

    // sets the corrected sequence separation
    //
    static void set_seq_sep(int i)
    {
        if (i >= cfg::params::seq_sep || instance()._forcing)
        {
            instance()._seq_sep = i;
        }

        else
        {
            instance()._seq_sep = cfg::params::seq_sep;
        }
    }

    // sets the network policy from direct value
    //
    static void set_net_policy(policy net_policy)
    {
        instance()._net_policy = net_policy;
    }

    // parses the network policy from a string; if error then throws an exception.
    //
    static void set_net_policy(std::string const& net_policy)
    {
        if (net_policy == "ca")
        {
            instance()._net_policy = policy::CA;
        }

        else if (net_policy == "cb")
        {
            instance()._net_policy = policy::CB;
        }

        else if (net_policy == "closest")
        {
            instance()._net_policy = policy::CLOSEST;
        }

        else
        {
            throw std::runtime_error("incorrect network policy argument: \"" + net_policy + "\"");
        }
    }

    // sets the bond control from direct value.
    //
    static void set_bond_control(control bond_control)
    {
        instance()._bond_control = bond_control;
    }

    // parses the bond control from a string; if error then throws an exception.
    //
    static void set_bond_control(std::string const& bond_control)
    {
        if (bond_control == "strict")
        {
            instance()._bond_control = control::STRICT_BOND;
        }

        else if (bond_control == "weak")
        {
            instance()._bond_control = control::WEAK_BOND;
        }

        else
        {
            throw std::runtime_error("incorrect bond control argument: \"" + bond_control + "\"");
        }
    }

    // sets the interaction type from direct value.
    //
    static void set_interaction_type(interaction_type interaction_type)
    {
        instance()._interaction_type = interaction_type;
    }

    // parses the interaction type from a string; if error then throws an exception.
    //
    static void set_interaction_type(std::string const& interaction_type)
    {
        parameters& params = instance();
        if (interaction_type == "all")
        {
            params._interaction_type = interaction_type::ALL;
        }

        else if (interaction_type == "multiple")
        {
            params._interaction_type = interaction_type::MULTIPLE;
        }

        else if (interaction_type == "one")
        {
            params._interaction_type = interaction_type::ONE;
        }

        else
        {
            throw std::runtime_error("incorret interaction type argument: \"" + interaction_type + "\"");
        }
    }

    // sets the forcing flag
    // FIXME this has effect only if it is called before any other distance setter... but it would be overkill to
    // have a builder pattern just because of this little one. proposed solution: just add this as an explicit ctor param.
    //
    static void set_forcing(bool flag)
    {
        instance()._forcing = flag;
    }

    static void set_hbond_realistic(bool flag)
    {
        instance()._hbond_realistic = flag;
    }

    // sets the path of the pdb file that has to be read
    //
    static void set_pdb_path(std::filesystem::path const& path)
    {
        instance()._pdb_path = path;
    }

    // sets the path of the generated rin
    //
    static void set_out_path(std::filesystem::path const& path)
    {
        instance()._output_path = path;
    }

    // sets the path of the logfile
    //
    static void set_log_path(std::filesystem::path const& path)
    {
        instance()._log_path = path;
    }

public:
    // the various getters
    //
    [[nodiscard]]
    static double get_distance_h()
    {
        return instance()._hbond_dist;
    }

    [[nodiscard]]
    static double get_distance_vdw()
    {
        return instance()._vdw_dist;
    }

    [[nodiscard]]
    static double get_distance_ionic()
    {
        return instance()._ionic_dist;
    }

    [[nodiscard]]
    static double get_distance_generic()
    {
        return instance()._generic_dist;
    }

    [[nodiscard]]
    static double get_distance_pication()
    {
        return instance()._pication_dist;
    }

    [[nodiscard]]
    static double get_distance_pipistack()
    {
        return instance()._pipistack_dist;
    }

    [[nodiscard]]
    static int get_seq_sep()
    {
        return instance()._seq_sep;
    }

    [[nodiscard]]
    static bool get_hbond_realistic()
    {
        return instance()._hbond_realistic;
    }

    [[nodiscard]]
    static policy get_net_policy()
    {
        return instance()._net_policy;
    }

    [[nodiscard]]
    static interaction_type get_interaction_type()
    {
        return instance()._interaction_type;
    }

    [[nodiscard]]
    static control get_bond_control()
    {
        return instance()._bond_control;
    }

    [[nodiscard]]
    static std::string get_pdb_name()
    {
        return instance()._pdb_path.stem().string();
    }

    [[nodiscard]]
    static std::filesystem::path const& get_pdb_path()
    {
        return instance()._pdb_path;
    }

    [[nodiscard]]
    static std::filesystem::path const& get_output_path()
    {
        return instance()._output_path;
    }

    [[nodiscard]]
    static std::filesystem::path const& get_log_path()
    {
        return instance()._log_path;
    }

public:
    // prettify the runtime configuration
    //
    [[nodiscard]]
    static std::string pretty();
};