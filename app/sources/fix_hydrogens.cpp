#include "fix_hydrogens.h"

#include "config.h"
#include "log_manager.h"

#include <iostream>
#include <filesystem>

using lm = log_manager;
namespace fs = std::filesystem;

void fix_hydrogens(gemmi::Structure& structure, gemmi::HydrogenChange what)
{
#ifndef _MSC_VER
    // todo: is this necessary?
    if (structure.models.empty() || structure.models[0].chains.empty())
        throw std::runtime_error{"the input file is empty"};

    lm::main()->info("fixing hydrogens...");

    gemmi::setup_entities(structure);
    auto h1 = gemmi::count_hydrogen_sites(structure);

    // todo add non-linux alternatives
    fs::path home{getenv("HOME")};

    if (what == gemmi::HydrogenChange::Remove)
        gemmi::remove_hydrogens(structure);
    else
    {
        auto const res_names = structure.models[0].get_all_residue_names();
        auto monlib = gemmi::read_monomer_lib(
            (home / std::filesystem::path{cfg::monomer_lib_dir}).string(),
            res_names,
            gemmi::cif::read_file,
            {},
            true
        );

        for (size_t i = 0; i < structure.models.size(); ++i)
            prepare_topology(structure, monlib, i, what, false, &std::cerr);
    }

    auto h2 = gemmi::count_hydrogen_sites(structure);

    lm::main()->info("hydrogen count before fix: {} after fix: {}", h1, h2);
#endif
}