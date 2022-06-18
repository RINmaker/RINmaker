#pragma once

#include <string>
#include <utility>
#include "ns_record.h"

namespace chemical_entity
{
class aminoacid;
}

namespace secondary_structure
{
class base
{
public:
    [[nodiscard]]
    virtual std::string pretty_with(chemical_entity::aminoacid const&) const;
};

class loop : public base
{
public:
    [[nodiscard]]
    std::string pretty_with(chemical_entity::aminoacid const&) const override;
};

class sheet_piece : public base
{
private:
    record::sheet_piece const _record;

public:
    explicit sheet_piece(record::sheet_piece const& record);

    [[nodiscard]]
    std::string pretty_with(chemical_entity::aminoacid const&) const override;
};

class helix : public base
{
private:
    record::helix const _record;

public:
    explicit helix(record::helix const& record);

    [[nodiscard]]
    std::string pretty_with(chemical_entity::aminoacid const&) const override;
};
}
