#include "ns_secondary_structure.h"

#include "ns_record.h"
#include "ns_chemical_entity.h"

using namespace secondary_structure;

sheet_piece::sheet_piece(record::sheet_piece const& record) :
        base(), _record(record)
{}

helix::helix(record::helix const& record) :
        base(), _record(record)
{}

std::string sheet_piece::pretty_with(chemical_entity::aminoacid const& res) const
{
    return "SHEET:" +
           _record.get_id() + ":" +
           std::to_string(1 + res.sequence_number() - _record.init_seq_number());
}

std::string helix::pretty_with(chemical_entity::aminoacid const& res) const
{
    return "HELIX:" +
           std::to_string(_record.serial()) + ":" +
           std::to_string(1 + res.sequence_number() - _record.init_seq_number());
}

std::string loop::pretty_with(chemical_entity::aminoacid const&) const
{ return "LOOP"; } // TODO config

std::string base::pretty_with(chemical_entity::aminoacid const&) const
{ return "NONE"; } // TODO config