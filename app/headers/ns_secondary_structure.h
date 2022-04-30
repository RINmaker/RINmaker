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
    record::sheet_piece const _record;
    std::string const _sheet_id, _seq_dist;

public:
    sheet_piece(record::sheet_piece record, chemical_entity::aminoacid const& res);

    [[nodiscard]]
    std::string pretty() const override;
};

class helix : public base
{
private:
    record::helix const _record;
    std::string const _helix_serial, _seq_dist;

public:
    helix(record::helix record, chemical_entity::aminoacid const& res);

    [[nodiscard]]
    std::string pretty() const override;
};
}
