#include "secondary_structures.h"

#include "ns_record.h"
#include "ns_chemical_entity.h"

using namespace structure;

std::string sheet_piece::pretty() const
{
    return "SHEET:" + _record.get_id() + ":" + std::to_string(1 + _res.sequence_number() - _record.init_seq_number());
}

std::string helix::pretty() const
{
    return "HELIX:" + std::to_string(_record.serial()) + ":" + std::to_string(1 + _res.sequence_number() - _record.init_seq_number());
}

std::string loop::pretty() const { return "LOOP"; } // TODO config

std::string base::pretty() const { return "NONE"; } // TODO config