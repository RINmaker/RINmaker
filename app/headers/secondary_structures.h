#pragma once

#include <string>
#include <utility>
#include "ns_record.h"

namespace chemical_entity
{
class aminoacid;
}

namespace structure
{
class base
{
public:
    [[nodiscard]]
    virtual std::string pretty() const;
};

class loop : public base
{
public:
    [[nodiscard]]
    std::string pretty() const override;
};

class sheet_piece : public base
{
private:
    records::sheet_piece const _record;
    chemical_entity::aminoacid const& _res;

public:
    sheet_piece(records::sheet_piece record, chemical_entity::aminoacid const& res) :
            base(),
            _record(std::move(record)),
            _res(res)
    {}

    [[nodiscard]]
    std::string pretty() const override;
};

class helix : public base
{
private:
    records::helix const _record;
    chemical_entity::aminoacid const& _res;

public:
    helix(records::helix record, chemical_entity::aminoacid const& res) :
            base(),
            _record(std::move(record)),
            _res(res)
    {}

    [[nodiscard]]
    std::string pretty() const override;
};
}
