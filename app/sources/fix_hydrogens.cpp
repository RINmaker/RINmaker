#include "fix_hydrogens.h"

#include "config.h"
#include "log_manager.h"

#include <optional>
#include <filesystem>

using lm = log_manager;
namespace fs = std::filesystem;

#ifdef _MSC_VER
#include <windows.h>
fs::path exec_path()
{
    TCHAR  buffer[MAX_PATH] = { 0 };
    GetModuleFileName(NULL, buffer, MAX_PATH);
    std::string path_str{buffer};
    return { path_str };
}
#else
fs::path exec_path()
{
    return fs::canonical({ "/proc/self/exe" });
}
#endif

class custom_stringbuf final : public std::stringbuf
{
public:
    int sync() override
    {
        // Get contents.
        auto s = this->str();

        // Compact in one line.
        replace_if(s.begin(), s.end(), [](auto ch) { return ch == '\n'; }, ' ');

        // Log...
        lm::main()->warn("[fix hydrogens] gemmi warning: {}", s);

        // Clear buffer.
        this->str("");
        return 0;
    }
};

void fix_hydrogens(gemmi::Structure& structure, gemmi::HydrogenChange what)
{
    if (structure.models.empty() || structure.models.front().chains.empty())
    {
        lm::main()->warn("[fix hydrogens] the input file has no models/has empty models!");
        return;
    }

    gemmi::setup_entities(structure);
    auto h1 = gemmi::count_hydrogen_sites(structure);

    if (what == gemmi::HydrogenChange::Remove)
    {
        lm::main()->info("[fix hydrogens] removing hydrogens...");
        gemmi::remove_hydrogens(structure);
    }
    else
    {
        auto const res_names = structure.models.front().get_all_residue_names();

        fs::path monomer_dir = exec_path().parent_path() / std::filesystem::path{ cfg::monomer_lib_name };

        lm::main()->info("[fix hydrogens] reading monomer library at: {}", monomer_dir.string());

        std::optional<gemmi::MonLib> monlib{std::nullopt};
        try
        {
            monlib = gemmi::read_monomer_lib(
                monomer_dir.string(),
                res_names,
                gemmi::cif::read_file,
                {},
                true
            );
        }
        catch (std::exception const& e)
        {
            lm::main()->error("[fix hydrogens] gemmi threw an exception: {}", e.what());
            return;
        }

        custom_stringbuf buff{};
        std::ostream warn{&buff};
        for (size_t i = 0; i < structure.models.size(); ++i)
        {
            lm::main()->info("[fix hydrogens] processing model {} of {}...", i+1, structure.models.size());
            prepare_topology(structure, *monlib, i, what, false, &warn);
        }

    }

    auto h2 = gemmi::count_hydrogen_sites(structure);
    lm::main()->info(
        "[fix hydrogens] hydrogen count (summing all {} models) before: {}, after: {}",
        structure.models.size(), h1, h2);
}