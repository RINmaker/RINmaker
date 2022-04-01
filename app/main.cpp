#include "utils.h"

using namespace std;
using lm = log_manager;

int main(int argc, const char* argv[]) {
    try {
        int i = 4;
        if (readArgs(argc, argv)) {

            lm::main()->debug("path to PDB input file: " + parameters::get_pdb_path().string());
            lm::main()->debug("path to output xml file: " + parameters::get_output_path().string());

            // TODO write to cout and also log
            // does this hold? --lore
            lm::main()->info("params summary: " + parameters::pretty());

            auto rm = rin_maker::rin_maker(parameters::get_pdb_path());
            rm(parameters::get_interaction_type(), parameters::get_net_policy()).consume_to_xml();

#           if _MSC_VER
            spdlog::drop_all();
#           endif
        }
    }

    catch (spdlog::spdlog_ex& e) {
        cerr << "log initialization failed: " << e.what() << endl;
        return 1;
    }

    catch (runtime_error& e) {
        cerr << "exception caught: " << e.what() << endl;
        return 1;
    }
}
