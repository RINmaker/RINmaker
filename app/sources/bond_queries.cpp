#include "bond_queries.h"

#include "noncovalent_bonds.h"
#include "rin_network.h"

#include "rin_params.h"

#include "chemical_entity.h"
#include "spatial/kdpoint.h"

using namespace bondfunctors;

void vdw::operator()(chemical_entity::atom const& a, chemical_entity::atom const& b)
{
    if (a.res().satisfies_minimum_separation(b.res()) &&
        a.distance(b) - (a.vdw_radius() + b.vdw_radius()) <= rin::parameters::global::instance().get().surface_dist_vdw())
    {
        string const unique_key = prelude::sort(a.res().id(), b.res().id());
        if (_bonded.find(unique_key) == _bonded.end())
        {
            // TODO verbosity: bond accepted
            _net.new_bond<bonds::vdw>(a, b);
            ++_nbonds;

            _bonded.insert(unique_key);
        }

        // TODO verbosity: bond rejected
    }
}

void ionic::operator()(chemical_entity::ionic_group const& a, chemical_entity::ionic_group const& b)
{
    if (a.res().satisfies_minimum_separation(b.res()) && a.charge() == -b.charge())
    {
        // TODO verbosity: bond accepted
        _net.new_bond<bonds::ionic>(a, b);
        ++_nbonds;
    }

    // TODO verbosity: bond rejected
}

void hydrogen::operator()(chemical_entity::atom const& acceptor, chemical_entity::atom const& donor)
{
    if (acceptor.res().satisfies_minimum_separation(donor.res()))
    {
        if (!(acceptor.res() == donor.res()))
        {
            auto hydrogens = donor.attached_hydrogens();
            for (auto* h: hydrogens)
            {
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


// FIXME possibile bug
// TODO mark resolved?
void pication::operator()(chemical_entity::atom const& cation, chemical_entity::ring const& ring)
{
    if (ring.res().satisfies_minimum_separation(cation.res(), rin::parameters::global::instance().get().sequence_separation()))
    {
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

void pipistack::operator()(chemical_entity::ring const& a, chemical_entity::ring const& b)
{
    double nc1 = a.angle_between_normal_and_centres_joining(b);
    double nc2 = b.angle_between_normal_and_centres_joining(a);
    double nn = a.angle_between_normals(b);
    double mn = a.closest_distance_between_atoms(b);

    if (a.res().satisfies_minimum_separation(b.res()) &&
        (0 <= nn && nn <= cfg::params::pipistack_normal_normal_angle_range) &&
        ((0 <= nc1 && nc1 <= cfg::params::pipistack_normal_centre_angle_range) ||
         (0 <= nc2 && nc2 <= cfg::params::pipistack_normal_centre_angle_range)) &&
        mn <= cfg::params::max_pipi_atom_atom_distance)
    {
        std::string const unique_key = prelude::sort(a.res().id(), b.res().id());
        if (_bonded.find(unique_key) == _bonded.end())
        {
            // TODO verbosity: bond accepted
            _net.new_bond<bonds::pipistack>(a, b, nn);
            ++_nbonds;

            _bonded.insert(unique_key);
        }

        // TODO verbosity: bond rejected
    }
}

void generico::operator()(chemical_entity::atom const& a, chemical_entity::atom const& b)
{
    if (a.res().satisfies_minimum_separation(b.res()))
    {
        string const unique_key = prelude::sort(a.res().id(), b.res().id());
        if (_bonded.find(unique_key) == _bonded.end())
        {
            // TODO verbosity: bond accepted
            _net.new_bond<bonds::generico>(a.res(), b.res());
            ++_nbonds;

            _bonded.insert(unique_key);
        }

        // TODO verbosity: bond rejected?
    }
}