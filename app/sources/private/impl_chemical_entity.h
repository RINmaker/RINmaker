#pragma once

#include "ns_chemical_entity.h"

#include <string>
#include <memory>
#include <vector>
#include <optional>
#include <array>

#include "ns_record.h"
#include "ns_secondary_structure.h"

struct chemical_entity::aminoacid::impl final
{
public:
    std::vector<chemical_entity::atom> _atoms;

    std::optional<chemical_entity::ring> primary_ring = std::nullopt, secondary_ring = std::nullopt;
    std::optional<chemical_entity::ionic_group> positive_ionic_group = std::nullopt, negative_ionic_group = std::nullopt;

    std::unique_ptr<secondary_structure::base> secondary_structure{std::make_unique<secondary_structure::base>()};

    std::optional<chemical_entity::atom> alpha_carbon = std::nullopt;
    std::optional<chemical_entity::atom> beta_carbon = std::nullopt;

    std::string chain_id;

    std::string name;

    int sequence_number = 0;

    std::string id;

    std::string pdb_name;

    std::array<double, 3> pos;
};

struct chemical_entity::aminoacid::component::impl final
{
public:
    std::weak_ptr<aminoacid::impl> res_impl;

    explicit impl(aminoacid const& res) : res_impl{res.pimpl}
    {}
};

struct chemical_entity::atom::impl final
{
public:
    record::atom record;
};

struct chemical_entity::ring::impl final
{
public:
    std::vector<chemical_entity::atom> atoms;

    std::array<double, 3> normal{};
    double mean_radius;
};

struct chemical_entity::ionic_group::impl
{
public:
    std::vector<chemical_entity::atom> const atoms;
    int const charge;
};