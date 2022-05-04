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
            lm::main()->debug("path to output xml file: " + parsed_args.out_dir.string());

            // TODO write to cout and also log
            // does this hold? --lore
            lm::main()->info("params summary: " + parsed_args.params.pretty());

            auto const models = rin::maker::parse_models(parsed_args.pdb_path);

            for (size_t i = 0; i < models.size(); ++i)
            {
                auto out_file = parsed_args.pdb_path.filename();
                out_file.replace_filename(out_file.stem().string() + "_" + to_string(i+1));
                out_file.replace_extension(".graphml");

                // create rin::maker, create graph and write to graphml
                models[i]()(parsed_args.params).write_to_file(parsed_args.out_dir / out_file);
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
