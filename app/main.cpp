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

        // this might change if reduce has been invoked
        auto final_input = parsed_args.input;
        if (parsed_args.reduce)
        {
            lm::console()->info("calling reduce...");

            if (parsed_args.input.extension() == ".pdb")
            {
                // file path for the reduced input file
                auto reduced_input =
                    parsed_args.input.parent_path() / std::filesystem::path{parsed_args.input.stem().string() + "_reduced.pdb"};

                // call reduce on the file with no console output; redirect output to the new file
                int success = std::system(("reduce -Quiet " + parsed_args.input.string() + " > " + reduced_input.string()).c_str());

                // check exit code
                if (success != 0)
                    lm::console()->error("calling reduce returned a nonzero exit code. the original input file will be used.");
                else
                {
                    lm::console()->info("calling reduce was successful");
                    final_input = reduced_input;
                }
            }
            else
                log_manager::console()->warn("input file is not a pdb, calling reduce aborted.");
        }

        lm::console()->info("final input file={}", final_input.filename().string());

        std::optional<gemmi::Structure> maybe_protein;

        // parse file
        if (final_input.extension() == ".pdb")
            maybe_protein = gemmi::read_pdb_file(final_input);
        else if (final_input.extension() == ".cif")
        {
            auto document = gemmi::cif::read_file(final_input);
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
