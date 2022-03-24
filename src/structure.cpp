#include "structure.h"

#include "records.h"
#include "entity.h"

using namespace structure;

std::string sheet_piece::pretty_with(entities::aminoacid const& res) const
{
    return "SHEET:" + _record.get_id() + ":" + std::to_string(1 + res.sequence_number() - _record.init_seq_number());
}

std::string helix::pretty_with(entities::aminoacid const& res) const
{
    return "HELIX:" + std::to_string(_record.serial()) + ":" + std::to_string(1 + res.sequence_number() - _record.init_seq_number());
}

std::string loop::pretty_with(entities::aminoacid const&) const { return "LOOP"; } // TODO config

std::string base::pretty_with(entities::aminoacid const&) const { return "NONE"; } // TODO config