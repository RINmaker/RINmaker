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
            std::filesystem::create_directory(parsed.log_path);
            log_manager::initialize(parsed.log_path);

            lm::main()->debug("path to PDB input file: " + parsed.pdb_path.string());
            lm::main()->debug("path to output xml file: " + parsed.out_path.string());

            // TODO write to cout and also log
            // does this hold? --lore
            // lm::main()->info("params summary: " + parameters::pretty());

            // parse file and build acceleration structures
            auto rm = rin::maker(parsed.pdb_path);
            rm(parsed.params).consume_to_xml();

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
