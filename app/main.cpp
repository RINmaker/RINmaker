#include "utils.h"
#include "rin_params.h"

using namespace std;
using lm = log_manager;
int main(int argc, const char* argv[])
{
    try
    {
        optional<arguments> maybe_args;
        if (read_args(argc, argv, maybe_args) && maybe_args.has_value())
        {
            arguments parsed = maybe_args.value();

            lm::main()->debug("path to PDB input file: " + parsed.pdb_path.string());
            lm::main()->debug("path to output xml file: " + parsed.out_path.string());

            // TODO write to cout and also log
            // does this hold? --lore
            lm::main()->info("params summary: " + parsed.params.pretty());

            // parse file and build acceleration structures
            auto rm = rin::maker(parsed.pdb_path);

            // create rin and write to graphml
            rm(parsed.params).write_to_file(parsed.params, parsed.out_path);

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
