#pragma warning(push, 0)

#include <gemmi/pdb.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>

#pragma warning(pop)

#include <cstdlib>

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
        if (!maybe_args.has_value())
            return 0;

        // the actual command line arguments
        auto const& parsed_args = *maybe_args;

        std::optional<gemmi::Structure> maybe_protein;

        // parse file
        if (parsed_args.input.extension() == ".pdb")
            maybe_protein = gemmi::read_pdb_file(parsed_args.input);
        else if (parsed_args.input.extension() == ".cif")
        {
            auto document = gemmi::cif::read_file(parsed_args.input);
            maybe_protein = gemmi::make_structure(document);
        }

        if (!maybe_protein.has_value())
        {
            lm::console()->error("only .pdb and .cif file formats are supported!");
            return 1;
        }

        auto protein = *maybe_protein;

        lm::console()->info("found {} models in protein {}", protein.models.size(), protein.name);

        // do all models
        if (holds_alternative<output_directory>(parsed_args.output))
        {
            lm::console()->info("selected all models, computing...");

            auto dir = get<output_directory>(parsed_args.output).value;
            for (auto const& model: protein.models)
            {
                std::filesystem::path file = protein.name + "_" + model.name + ".graphml";

                // create rin::maker, create graph and write to graphml
                rin::maker{model, protein, parsed_args.skip_water}(parsed_args.params).write_to_file(dir / file);
            }
        }

        // otherwise do just the first one
        else if (holds_alternative<output_file>(parsed_args.output))
        {
            lm::console()->info("selected first model, computing...");
            auto const file = get<output_file>(parsed_args.output).value;

            // create rin::maker, create graph and write to graphml
            rin::maker{protein.first_model(), protein, parsed_args.skip_water}(parsed_args.params).write_to_file(file);
        }

        lm::console()->info("done");

#       if _MSC_VER
        spdlog::drop_all();
#       endif
    }

    catch (spdlog::spdlog_ex const& e)
    {
        cerr << "log initialization failed: " << e.what() << endl;
        return 1;
    }

    catch (exception const& e)
    {
        cerr << "exception caught: " << e.what() << endl;
        return 1;
    }
}
