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

            //lm::main()->debug("path to PDB input file: " + parsed_args.pdb_path.string());
            //lm::main()->debug("path to output xml file: " + parsed_args.out_dir.string());

            // TODO write to cout and also log
            // does this hold? --lore
            //lm::main()->info("params summary: " + parsed_args.params.pretty());

            //lm::console()->info("file={}", parsed_args.input.filename().string());

            // parse file
            auto protein = gemmi::read_pdb_file(parsed_args.input);

            //lm::console()->info("found {} models in protein {}", protein_structure.models.size(), protein_structure.name);

            // do all models
            if (holds_alternative<output_directory>(parsed_args.output))
            {
                auto dir = get<output_directory>(parsed_args.output).value;
                for (auto const& model: protein.models)
                {
                    std::filesystem::path file = protein.name + "_" + model.name + ".graphml";

                    // create rin::maker, create graph and write to graphml
                    rin::maker{model, protein}(parsed_args.params).write_to_file(dir / file);
                }
            }

            // otherwise do just the first one
            else if (holds_alternative<output_file>(parsed_args.output))
            {
                auto const file = get<output_file>(parsed_args.output).value;
                // create rin::maker, create graph and write to graphml
                rin::maker{protein.first_model(), protein}(parsed_args.params).write_to_file(file);
            }
            //lm::console()->info("done");

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
