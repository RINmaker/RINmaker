#pragma once

#include "ns_chemical_entity.h"

#include <string>
#include <memory>
#include <vector>
#include <optional>
#include <array>

struct chemical_entity::aminoacid::impl final
{
public:
    std::vector<chemical_entity::atom> atoms;

    std::optional<chemical_entity::ring> primary_ring = std::nullopt;
    std::optional<chemical_entity::ring> secondary_ring = std::nullopt;
    std::optional<chemical_entity::ionic_group> positive_ionic_group = std::nullopt;
    std::optional<chemical_entity::ionic_group> negative_ionic_group = std::nullopt;

    std::string secondary_structure_name;

    std::optional<chemical_entity::atom> alpha_carbon = std::nullopt;
    std::optional<chemical_entity::atom> beta_carbon = std::nullopt;

    std::string chain_id;

    std::string name;

    int sequence_number = 0;

    std::string id;

    std::string protein_name;

    std::array<double, 3> position;
};

struct chemical_entity::atom::impl final
{
public:
    gemmi::Atom record;

    explicit impl(gemmi::Atom  record) : record{std::move(record)}
    {}
};

struct chemical_entity::ring::impl final
{
public:
    std::vector<chemical_entity::atom> atoms;

    std::array<double, 3> normal;
};

struct chemical_entity::ionic_group::impl
{
public:
    std::vector<chemical_entity::atom> const atoms;
    int const charge;

    impl(std::vector<chemical_entity::atom>  atoms, int const charge) : atoms{std::move(atoms)}, charge{charge}
    {}
};