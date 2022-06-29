#pragma warning(push, 0)

#include <gemmi/pdb.hpp>

#pragma warning(pop)

#include "cli_utils.h"

#include "rin_maker.h"
#include "rin_params.h"

using namespace std;
using lm = log_manager;
int main(int argc, const char* argv[])
{
    try
    {
        optional<arguments> maybe_args = read_args(argc, argv);
        if (maybe_args.has_value())
        {
            auto const parsed_args = *maybe_args;

            lm::main()->debug("path to PDB input file: " + parsed_args.pdb_path.string());
            lm::main()->debug("path to output xml file: " + parsed_args.out_dir.string());

            // TODO write to cout and also log
            // does this hold? --lore
            lm::main()->info("params summary: " + parsed_args.params.pretty());

            lm::console()->info("file={}", parsed_args.pdb_path.filename().string());

            auto protein_structure = gemmi::read_pdb_file(parsed_args.pdb_path);

            lm::console()->info("found {} models in protein {}", protein_structure.models.size(), protein_structure.name);

            for (auto const& model: protein_structure.models)
            {
                std::filesystem::path out_file = protein_structure.name + "_" + model.name + ".graphml";

                // create rin::maker, create graph and write to graphml
                rin::maker{model, protein_structure}(parsed_args.params).write_to_file(parsed_args.out_dir / out_file);
                lm::console()->info("model={}", model.name);
            }
            lm::console()->info("done");

#           if _MSC_VER
            spdlog::drop_all();
#           endif
        }
    }

    catch (spdlog::spdlog_ex& e)
    {
        cerr << "log initialization failed: " << e.what() << endl;
        return 1;
    }

    catch (exception& e)
    {
        cerr << "exception caught: " << e.what() << endl;
        return 1;
    }
}
