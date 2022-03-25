#include "utils.h"

using namespace std;
using lm = log_manager;

int main(int argc, const char* argv[]) {
    try {
        if (readArgs(argc, argv)) {

            lm::main()->debug("path to PDB input file: " + parameters::get_pdb_path().string());
            lm::main()->debug("path to output xml file: " + parameters::get_output_path().string());

            // TODO write to cout and also log
            // does this hold? --lore
            lm::main()->info("params summary: " + parameters::pretty());

            auto rm = rin_maker::build(parameters::get_net_policy(), parameters::get_pdb_path());
            if (rm == nullptr)
                throw std::runtime_error("unknown bonds policy");

            rm->get_graph().consume_to_xml();
            delete rm;

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
