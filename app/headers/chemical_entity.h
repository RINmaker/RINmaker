#pragma once

#include <string>
#include <vector>
#include <array>
#include <tuple>

#include "config.h"

#include "runtime_params.h"

#include "graphml_output.h"

#include "energy.h"

#include "spatial/geometry.h"
#include "prelude.h"
#include "spatial/kdpoint.h"

#include "pdb_records.h"
#include "secondary_structures.h"

#include "log_manager.h"

class pdb_data;
namespace rin_maker {class base; struct rin_maker; }

namespace chemical_entity {

    class atom;
    class ring;
    class ionic_group;

    class aminoacid : public kdpoint<3> {
    private:
        std::vector<atom const *> _atoms = std::vector<atom const *>();
        atom const *_alpha_carbon = nullptr;
        atom const *_beta_carbon = nullptr;

        ring const *_primary_ring = nullptr;
        ring const *_secondary_ring = nullptr;

        ionic_group const *_positive_ionic_group = nullptr;
        ionic_group const *_negative_ionic_group = nullptr;

    private:
        /**
         * The ID of the aminoacid's <b>chain</b>.
         * <br>
         * Examples: <i>A</i>, <i>B</i>, <i>C</i>
         */
        std::string _chain_id;

        /**
         * The name of the aminoacid.
         * <br>
         * Examples: <i>PHE</i>, <i>HYS</i>, <i>TRP</i>
         */
        std::string _name;

        /**
         * The sequence number of the aminoacid.
         */
        int _sequence_number = 0;

        /**
         * A <b>unique</b> identifier for the aminoacid.
         * <br>
         * Format (from Padua): <i>CHAIN:SEQ:_:NAME</i>
         */
        std::string _id;

        /**
         * The secondary structure can be either <i>LOOP</i>, <i>HELIX</i> or <i>SHEET</i>
         * when information is available in the PDB, otherwise <i>NONE</i>.
         */
        structure::base *_secondary_structure = new structure::base();

    private:
        friend class ::pdb_data;
        friend class rin_maker::base;
        friend struct rin_maker::rin_maker;

        // lifetime is managed by pdb_data
        explicit aminoacid(std::vector<records::atom> const &records);

        ~aminoacid();

    public:
        /**
         * @return A collection of all the aminoacid's atoms.
         */
        [[nodiscard]]
        std::vector<atom const *> const &atoms() const { return _atoms; }

        /**
         * @return A pointer to the <b>alpha</b> backbone.
         */
        [[nodiscard]]
        atom const *ca() const { return _alpha_carbon; }

        /**
         * @return A pointer to the <b>beta</b> backbone.
         */
        [[nodiscard]]
        atom const *cb() const { return _beta_carbon; }

        /**
         * @return A pointer to the <b>primary</b> aromatic ring.
         */
        [[nodiscard]]
        ring const *primary_ring() const { return _primary_ring; }

        /**
         * @return A pointer to the <b>secondary</b> aromatic ring.
         */
        [[nodiscard]]
        ring const *secondary_ring() const { return _secondary_ring; }

        /**
         * @return A pointer to the <b>positive</b> ionic group.
         */
        [[nodiscard]]
        ionic_group const *positive_ionic_group() const { return _positive_ionic_group; }

        /**
         * @return A pointer to the <b>negative</b> ionic group.
         */
        [[nodiscard]]
        ionic_group const *negative_ionic_group() const { return _negative_ionic_group; }

    public:
        /**
         * @return The name of the aminoacid.
         *
         * Examples: <i>PHE</i>, <i>HYS</i>, <i>TRP</i>
         */
        [[nodiscard]]
        std::string const &name() const { return _name; }

        /**
         * @return The ID of the aminoacid's <b>chain</b>.
         *
         * Examples: <i>A</i>, <i>B</i>, <i>C</i>
         */
        [[nodiscard]]
        std::string const &chain_id() const { return _chain_id; }

        /**
         * @return A <b>unique</b> identifier for the aminoacid.
         *
         * Format (from Padua): <i>CHAIN:SEQ:_:NAME</i>
         */
        [[nodiscard]]
        std::string const &id() const { return _id; }

        /**
         * @return The sequence number of the aminoacid.
         */
        [[nodiscard]]
        int sequence_number() const { return _sequence_number; }

    public:
        /**
         * @param rhs - Another aminoacid
         * @return True iff <i>this</i> and <i>rhs</i> have the same ID
         */
        bool operator==(aminoacid const &rhs) const { return _id == rhs._id; }

        /**
         * @param rhs Another aminoacid
         * @return True iff they do <b>not</b> have the same ID
         */
        bool operator!=(aminoacid const &rhs) const { return !(*this == rhs); }

        /**
         * @param aa Another aminoacid
         * @param minimum_separation Minimum accepted sequence separation (inclusive)
         * @return True if they are separated by (at least) <i>seq_sep</i>
         */
        [[nodiscard]]
        bool
        satisfies_minimum_separation(aminoacid const &aa, int minimum_separation = parameters::get_seq_sep()) const {
            if (*this == aa) {
                return false;
            }

            if (_chain_id != aa._chain_id) {
                return true;
            }

            return abs(_sequence_number - aa._sequence_number) >= minimum_separation;
        }

        // return a rin::node from this
        [[nodiscard]]
        rin::node to_node() const { return rin::node(*this); } // TODO pu� diventare un operatore di conversione

    public:
        // makes structure become <i>LOOP</i>
        void make_secondary_structure() {
            delete _secondary_structure;
            _secondary_structure = new structure::loop();
        }

        // makes structure become <i>HELIX</i>
        void make_secondary_structure(records::helix const &record) {
            delete _secondary_structure;
            _secondary_structure = new structure::helix(record);
        }

        // makes structure become <i>SHEET</i>
        void make_secondary_structure(records::sheet_piece const &record) {
            delete _secondary_structure;
            _secondary_structure = new structure::sheet_piece(record);
        }

        // ID of the secondary structure
        [[nodiscard]]
        std::string secondary_structure_id() const {
            return _secondary_structure->pretty_with(*this);
        }
    };

/**
 * Base class for the units that compose an aminoacid:
 * <br>
 * <i>atom</i>, <i>ring</i>, <i>ionic_group</i>
 */
    class component {
    protected:
        friend class aminoacid;

        // the aminoacids that "owns" this component.
        aminoacid const *_res;

        // we need to set its aminoacid
        explicit component(aminoacid const &res) : _res(&res) {}

        virtual ~component() = default;

    public:
        /**
         * @return A <b>reference</b> to its aminoacid.
         */
        [[nodiscard]]
        aminoacid const &res() const { return *_res; }
    };

/**
 * Represents an atom as both:
 * <ul>
 * <li>a component of aminoacid</li>
 * <li>a point in 3D space</li>
 * </ul>
 */
    class atom final : public kdpoint<3>, public component {
    private:
        // lifetime is managed by aminoacid
        friend class aminoacid;

        atom(records::atom const &record, aminoacid const &res)
                : kdpoint<3>({record.x(), record.y(), record.z()}), component(res), _record(record) {}

        ~atom() override = default;

    private:
        // the record it was parsed from
        records::atom _record;

    public:
        /**
         * @return The name of the atom.
         */
        [[nodiscard]]
        std::string const &name() const { return _record.name(); }

        /**
         * @return The chemical symbol of the atom.
         */
        [[nodiscard]]
        std::string const &symbol() const { return _record.element_name(); }

        /**
         * @return The temperature factor.
         */
        [[nodiscard]]
        double temp_factor() const { return _record.temp_factor(); }

        /**
         * @return The electric charge.
         */
        [[nodiscard]]
        int charge() const { return _record.charge(); }

        /**
         * @return True iff this atom <b>is in</b> a positive ionic group.
         */
        [[nodiscard]]
        bool is_in_a_positive_ionic_group() const;

        /**
         * @return True iff this atom <b>is in</b> a negative ionic group.
         */
        [[nodiscard]]
        bool is_in_a_negative_ionic_group() const;

        /**
         * @return True iff this atom is a donor of hydrogens.
         */
        [[nodiscard]]
        bool is_a_hydrogen_donor() const;

        int how_many_hydrogen_can_donate() const;

        /**
         * @return True iff this atom is an acceptor of hydrogens.
         */
        [[nodiscard]]
        bool is_a_hydrogen_acceptor() const;

        int how_many_hydrogen_can_accept() const;

        /**
         * @return A collection of its attached hydrogens. The collection might be empty.
         */
        [[nodiscard]]
        std::vector<atom const *> attached_hydrogens() const;

        /**
         * @return True if this atom <b>might</b> participate in a vdw bond.
         */
        [[nodiscard]]
        bool is_a_vdw_candidate() const;

        /**
         * @return Van der Waals radius of the atom.
         */
        [[nodiscard]]
        double vdw_radius() const;

        /**
         * @return The mass of the atom.
         */
        [[nodiscard]]
        double mass() const;

        /**
         * @return True iff this atom is a hydrogen.
         */
        [[nodiscard]]
        bool is_a_hydrogen() const { return symbol() == "H"; }

        /**
         * @return True iff atom is cation
         */
        [[nodiscard]]
        bool is_a_cation() const;

        /**
         * @return True iff atom is on main chain
         */
        [[nodiscard]]
        bool is_main_chain() const {
            return
                    this->name() == "C" ||
                    this->name() == "O" ||
                    this->name() == "H" ||
                    this->name() == "HA" ||
                    this->name() == "N";
        }
    };

/**
 * Represents an aromatic ring.
 */
    class ring final : public kdpoint<3>, public component {
    private:
        // atoms of the ring
        std::vector<atom const *> _atoms;

        // normal vector to the ring's plane
        std::array<double, 3> _normal{};

        // mean distance atom-centre
        double _mean_radius;

    private:
        friend class aminoacid;

        // only an aminoacid manages a ring's lifetime
        ring(std::vector<atom const *> const &atoms, aminoacid const &res);

        ~ring() override = default;

    public:
        /**
         * @return The normal to the plane of the ring.
         */
        [[nodiscard]]
        std::array<double, 3> const &normal() const { return _normal; }

        /**
         * @return The mean radius of the ring.
         */
        [[nodiscard]]
        double radius() const { return _mean_radius; }

        /**
         * @return Whether it could participate in a cation-pi bond.
         */
        [[nodiscard]]
        bool is_a_pication_candidate() const {
            string name = res().name();
            return name == "PHE" || name == "TYR" || (name == "TRP" && _atoms.size() == 6);
        }

        /**
         * @param other Another ring.
         * @return The distance between the closest atoms of the two rings.
         */
        [[nodiscard]]
        double closest_distance_between_atoms(ring const &other) const {
            double minimum = _atoms[0]->distance(*other._atoms[0]);
            for (auto *atom_1: _atoms) {
                for (auto *atom_2: other._atoms) {
                    double current = atom_1->distance(*atom_2);
                    if (current < minimum) {
                        minimum = current;
                    }
                }
            }

            return minimum;
        }

        /**
         * @param other Another ring.
         * @return The angle between the normals of the two rings, comprised between 0� and 90�.
         */
        [[nodiscard]]
        double angle_between_normals(ring const &other) const { return geom::d_angle<3>(_normal, other._normal); }

        /**
         * @param other Another ring.
         * @return The angle between this ring's normal and the centre-centre segment, comprised between 0� and 90�.
         */
        [[nodiscard]]
        double angle_between_normal_and_centres_joining(ring const &other) const {
            std::array<double, 3> const centres_joining((std::array<double, 3>) (*this - other));
            return geom::d_angle<3>(_normal, centres_joining);
        }

        /**
         * @param atom An atom.
         * @return An atom of this ring, the closest to other.
         */
        [[nodiscard]]
        atom const &atom_closest_to(atom const &atom) const;

        /**
         * @return The string representation of the ring.
         */
        [[nodiscard]]
        string name() const;
    };

/**
 * Represents a ionic group of an aminoacid.
 */
    class ionic_group final : public kdpoint<3>, public component {
    private:
        // all the atoms of the ionic group
        std::vector<atom const *> const _atoms;

        // its electric charge
        int const _charge;

    private:
        friend class aminoacid;

        // lifetime is managed by aminoacid
        ionic_group(std::vector<atom const *> const &atoms, int const &charge, aminoacid const &res);

        ~ionic_group() override = default;

    public:

        /**
         * @return The charge of the ionic group.
         */
        [[nodiscard]]
        int charge() const { return _charge; }

        /**
         * @return The ionion energy of the ionic group.
         */
        [[nodiscard]]
        double ionion_energy_q() const;

        /**
         * @return The string representation of the group.
         */
        [[nodiscard]]
        string name() const;
    };
}
