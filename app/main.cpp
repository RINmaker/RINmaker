#include <iostream>
#include <string>
#include <memory>
#include <cstdio>

#include <CLI/CLI.hpp>

#include "rin_maker.h"
#include "main.h"

#include "pdb_data.h"
#include "log_manager.h"

#include "config.h"
#include "runtime_params.h"

using namespace std;

int main(int argc, const char* argv[]) {
    try {
        if (prelude::readArgs(argc, argv)) {
            auto logger = log_manager::main();

            logger->debug("path to PDB input file: " + parameters::get_pdb_path().string());
            logger->debug("path to output xml file: " + parameters::get_output_path().string());
            logger->info("params summary: " + parameters::pretty()); // TODO scrivere in cout E ANCHE in log

            // TODO: remove use of smart ptr, useless.
            // unique_ptr<pdb_data> data = std::make_unique<pdb_data>();

            auto run = rin_maker::build(parameters::get_net_policy(), parameters::get_pdb_path());
            if (run == nullptr)
                throw std::runtime_error("unknown policy");
            run->get_graph().consume_to_xml();
            delete run;

#			if _MSC_VER
            spdlog::drop_all();
#			endif
        }
    }

    catch (spdlog::spdlog_ex& e) {
        cerr << "log initialization failed: " << e.what() << endl;
        return 1;
    }

    catch (runtime_error& e) {
        cerr << "exception caught: " << e.what() << ";; check logs" << endl;
        return 1;
    }
}
