#include "ns_secondary_structure.h"

#include "ns_record.h"
#include "ns_chemical_entity.h"

using namespace secondary_structure;

sheet_piece::sheet_piece(record::sheet_piece const& record, chemical_entity::aminoacid const& res) :
    base(),
    _record(record),
    _sheet_id{record.get_id()},
    _seq_dist{std::to_string(1 + res.sequence_number() - record.init_seq_number())}
{}

helix::helix(record::helix const& record, chemical_entity::aminoacid const& res) :
    base(),
    _record(record),
    _helix_serial{std::to_string(record.serial())},
    _seq_dist{std::to_string(1 + res.sequence_number() - record.init_seq_number())}
{}

std::string sheet_piece::pretty() const
{
    return "SHEET:" + _sheet_id + ":" + _seq_dist;
}

std::string helix::pretty() const
{
    return "HELIX:" + _helix_serial + ":" + _seq_dist;
}

std::string loop::pretty() const { return "LOOP"; } // TODO config

std::string base::pretty() const { return "NONE"; } // TODO config