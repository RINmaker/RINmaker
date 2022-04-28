#include "utils.h"

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
            auto const parsed_args = maybe_args.value();

            lm::main()->debug("path to PDB input file: " + parsed_args.pdb_path.string());
            lm::main()->debug("path to output xml file: " + parsed_args.out_path.string());

            // TODO write to cout and also log
            // does this hold? --lore
            lm::main()->info("params summary: " + parsed_args.params.pretty());

            auto const models = rin::maker::parse_models(parsed_args.pdb_path);

            int i = 0;
            for (auto const& rm : models)
            {
                // create rin and write to graphml
                auto const view = rm()(parsed_args.params);

                auto const out_path = parsed_args.out_path.parent_path() / (
                        parsed_args.pdb_path.stem().string() +
                        "_" +
                        to_string(i++) +
                        parsed_args.out_path.extension().string());
                view.write_to_file(out_path);
            }

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

    catch (runtime_error& e)
    {
        cerr << "exception caught: " << e.what() << endl;
        return 1;
    }
}
