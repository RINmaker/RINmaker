#pragma once

#include <array>
#include <string>
#include <set>
#include <tuple>

#include "noncovalent_bonds.h"
#include "chemical_entity.h"

#include "rin_network.h"
#include "log_manager.h"
#include "spatial/kdpoint.h"
#include "rin_params.h"

namespace bondfunctors {

class base {
protected:
    network& _net;
    int _nbonds;

protected:
    explicit base(network& net)
            : _net(net), _nbonds(0) {}

public:
    virtual ~base() = default;
};

/**
 * Van der Waals bond test.
 */
class vdw : public base {
private:
    // just to keep track of already accepted pairs; vdw is tested for ALL VDWS against ALL VDWS
    std::set<std::string> _bonded;

public:
    explicit vdw(network& net)
            : base(net) {}

    ~vdw() override { log_manager::main()->info("found {} vdw bonds", _nbonds); }

public:
    void operator()(chemical_entity::atom const& a, chemical_entity::atom const& b) {
        if (a.res().satisfies_minimum_separation(b.res()) &&
            a.distance(b) - (a.vdw_radius() + b.vdw_radius()) <= rin::parameters::global::instance().get().surface_dist_vdw()) {
            std::string const unique_key = prelude::sort(a.res().id(), b.res().id());
            if (_bonded.find(unique_key) == _bonded.end()) {
                // TODO verbosity: bond accepted
                _net.new_bond<bonds::vdw>(a, b);
                ++_nbonds;

                _bonded.insert(unique_key);
            }

            // TODO verbosity: bond rejected
        }
    }
};

/**
 * Ionic bond test.
 */
class ionic : public base {
public:
    explicit ionic(network& net)
            : base(net) {}

    ~ionic() override { log_manager::main()->info("found {} ionic bonds", _nbonds); }

public:
    // aminoacids must be separated and have opposite charges
    void operator()(chemical_entity::ionic_group const& a, chemical_entity::ionic_group const& b) {
        if (a.res().satisfies_minimum_separation(b.res()) && a.charge() == -b.charge()) {
            // TODO verbosity: bond accepted
            _net.new_bond<bonds::ionic>(a, b);
            ++_nbonds;
        }

        // TODO verbosity: bond rejected
    }
};

/**
 * Hydrogen bond test
 */
class hydrogen : public base {
public:
    explicit hydrogen(network& net)
            : base(net) {}

    ~hydrogen() override { log_manager::main()->info("found {} hydrogen bonds", _nbonds); }

public:
    // a bit complex, please refer to the paper
    void operator()(chemical_entity::atom const& acceptor, chemical_entity::atom const& donor) {
        if (acceptor.res().satisfies_minimum_separation(donor.res())) {
            if (!(acceptor.res() == donor.res())) {
                auto hydrogens = donor.attached_hydrogens();
                for (auto* h: hydrogens) {
                    std::array<double, 3> const da = (std::array<double, 3>) (acceptor - donor);
                    std::array<double, 3> const dh = (std::array<double, 3>) (*h - donor);
                    double angle_adh = geom::angle<3>(da, dh);

                    std::array<double, 3> const ha = (std::array<double, 3>) (acceptor - *h);
                    std::array<double, 3> const hd = (std::array<double, 3>) (donor - *h);
                    double angle_ahd = geom::angle<3>(ha, hd);

                    if (angle_adh <= cfg::params::hbond_angle) // 63
                    {
                        // TODO verbosity: bond accepted
                        _net.new_bond<bonds::hydrogen>(acceptor, donor, h, angle_ahd);
                        ++_nbonds;
                    }

                    // TODO verbosity: bond rejected
                }
            }
        }
    }
};

/**
 * PiCation bond test.
 */
class pication : public base {
public:
    explicit pication(network& net)
            : base(net) {}

    ~pication() override { log_manager::main()->info("found {} pication bonds", _nbonds); }

public:
    // FIXME possibile bug
    // TODO mark resolved?
    void operator()(chemical_entity::atom const& cation, chemical_entity::ring const& ring) {
        if (ring.res().satisfies_minimum_separation(cation.res(), rin::parameters::global::instance().get().sequence_separation())) {
            // TODO verbosity: bond accepted
            double theta = 90 - geom::d_angle<3>(ring.normal(), (std::array<double, 3>) (ring - cation));
            if (theta >= cfg::params::pication_angle) // 45
            {
                _net.new_bond<bonds::pication>(ring, cation, theta);
                ++_nbonds;
            }

            // TODO verbosity: bond rejected
        }
    }
};

/**
 * PiPiStacking bond test.
 */
class pipistack : public base {
private:
    // just to keep track of already accepted pairs; pipistack is tested for ALL RINGS against ALL RINGS
    std::set<std::string> _bonded;

public:
    explicit pipistack(network& net)
            : base(net) {}

    ~pipistack() override { log_manager::main()->info("found {} pipistack bonds", _nbonds); }

public:
    void operator()(chemical_entity::ring const& a, chemical_entity::ring const& b) {
        double nc1 = a.angle_between_normal_and_centres_joining(b);
        double nc2 = b.angle_between_normal_and_centres_joining(a);
        double nn = a.angle_between_normals(b);
        double mn = a.closest_distance_between_atoms(b);

        if (a.res().satisfies_minimum_separation(b.res()) &&
            (0 <= nn && nn <= cfg::params::pipistack_normal_normal_angle_range) &&
            ((0 <= nc1 && nc1 <= cfg::params::pipistack_normal_centre_angle_range) ||
             (0 <= nc2 && nc2 <= cfg::params::pipistack_normal_centre_angle_range)) &&
            mn <= cfg::params::max_pipi_atom_atom_distance) {
            std::string const unique_key = prelude::sort(a.res().id(), b.res().id());
            if (_bonded.find(unique_key) == _bonded.end()) {
                // TODO verbosity: bond accepted
                _net.new_bond<bonds::pipistack>(a, b, nn);
                ++_nbonds;

                _bonded.insert(unique_key);
            }

            // TODO verbosity: bond rejected
        }
    }
};

/**
 * Generic bond test (alpha-alpha or beta-beta).
 */
class generico : public base {
private:
    // just to keep track of already accepted pairs; generic is tested for ALL CARBONS against ALL CARBONS
    std::set<std::string> _bonded;

public:
    explicit generico(network& net)
            : base(net) {}

    ~generico() override { log_manager::main()->info("found {} generic bonds", _nbonds); }

public:
    // tecnicamente il test è sempre valido, ricordiamoci che questo è il test applicato ai vicini della range search!
    //
    void operator()(chemical_entity::atom const& a, chemical_entity::atom const& b) {
        if (a.res().satisfies_minimum_separation(b.res())) {
            string const unique_key = prelude::sort(a.res().id(), b.res().id());
            if (_bonded.find(unique_key) == _bonded.end()) {
                // TODO verbosity: bond accepted
                _net.new_bond<bonds::generico>(a.res(), b.res());
                ++_nbonds;

                _bonded.insert(unique_key);
            }
        }
    }
};
}
