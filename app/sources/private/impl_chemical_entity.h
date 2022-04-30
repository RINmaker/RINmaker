#pragma once

#include "ns_chemical_entity.h"

#include <string>
#include <memory>
#include <vector>
#include <array>

#include "ns_record.h"
#include "secondary_structures.h"

struct chemical_entity::aminoacid::impl final
{
public:
    std::vector<std::unique_ptr<chemical_entity::atom const>> _atoms;

    std::unique_ptr<chemical_entity::ring const> primary_ring, secondary_ring;
    std::unique_ptr<chemical_entity::ionic_group const> positive_ionic_group, negative_ionic_group;

    std::unique_ptr<structure::base> secondary_structure{std::make_unique<structure::base>()};

    chemical_entity::atom const* alpha_carbon = nullptr;
    chemical_entity::atom const* beta_carbon = nullptr;

    std::string chain_id;

    std::string name;

    int sequence_number = 0;

    std::string id;

    std::string pdb_name;
};

struct chemical_entity::aminoacid::component::impl final
{
public:
    std::weak_ptr<aminoacid::impl> res_impl;
    std::array<double, 3> res_pos;

    explicit impl(aminoacid const& res) : res_impl{res.pimpl}, res_pos{(std::array<double, 3>) res}
    {}
};

struct chemical_entity::atom::impl final
{
public:
    records::atom record;
};

struct chemical_entity::ring::impl final
{
public:
    std::vector<chemical_entity::atom const*> atoms;

    std::array<double, 3> normal{};
    double mean_radius;
};

struct chemical_entity::ionic_group::impl
{
public:
    std::vector<chemical_entity::atom const*> const atoms;
    int const charge;
};