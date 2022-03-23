#pragma once

#include <string>
#include "records.h"

// entity has to #include "structure.h"
// #include "entity.h"

namespace entities { class aminoacid; }

namespace structure
{
// 'base' secondary structure is the default for aminoacids
// it means 'no information'
//
class base
{
protected:
    friend class entities::aminoacid;

    base() = default;
    virtual ~base() = default;

public:
    // ritorna una stringa con le informazioni sulla struttura secondaria e l'amminoacido all'interno di essa
    //
    virtual std::string pretty_with(entities::aminoacid const&) const;
};

// 'loop' secondary structure
//
class loop : public base
{
private:
    friend class entities::aminoacid;

    loop() = default;
    ~loop() = default;

public:
    std::string pretty_with(entities::aminoacid const&) const;
};

// 'sheet piece' secondary structure
//
class sheet_piece : public base
{
private:
    friend class entities::aminoacid;
    records::sheet_piece const _record;

    sheet_piece(records::sheet_piece const& record) : base(), _record(record) {}
    ~sheet_piece() = default;

public:
    std::string pretty_with(entities::aminoacid const& res) const;
};

// 'helix' is secondary structure
//
class helix : public base
{
private:
    friend class entities::aminoacid;
    records::helix const _record;

    helix(records::helix const& record) : base(), _record(record) {}
    ~helix() = default;

public:
    std::string pretty_with(entities::aminoacid const& res) const;
};
}
