#pragma once

#include <vector>
#include <string>

#include <memory>
#include <functional>
#include <filesystem>

#include "rin_graph.h"
#include "ns_record.h"

namespace rin
{
struct parameters;

struct maker final
{
private:
    struct impl;
    std::shared_ptr<impl const> pimpl;

public:
    static std::vector<std::function<rin::maker(void)>> parse_models(std::filesystem::path const& pdb_path);

    maker(std::string const& pdb_name,
          std::vector<record::atom> const& atom_records,
          std::vector<record::ss> const& ssbond_records,
          std::vector<record::helix> const& helix_records,
          std::vector<record::sheet_piece> const& sheet_records);

    ~maker();

    rin::graph operator()(parameters const& params) const;
};
}