#pragma once

#include <vector>
#include <string>

#include <memory>
#include <filesystem>

#include "rin_graph.h"
#include "pdb_records.h"

namespace rin
{
struct parameters;

struct maker final
{
private:
    struct impl;
    std::shared_ptr<impl const> pimpl;

public:
    static std::vector<std::shared_ptr<rin::maker>> parse_models(std::filesystem::path const& pdb_path);

    maker(std::string const& pdb_name,
          std::vector<records::atom> const& atom_records,
          std::vector<records::ss> const& ssbond_records,
          std::vector<records::helix> const& helix_records,
          std::vector<records::sheet_piece> const& sheet_records);

    ~maker();

    rin::graph operator()(parameters const& params) const;
};
}