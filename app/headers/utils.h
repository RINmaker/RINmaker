#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <optional>
#include <filesystem>

#include <spdlog/spdlog.h>
#include <CLI/CLI.hpp>

#include "rin_params.h"
#include "log_manager.h"

struct arguments final
{
    rin::parameters params;
    std::filesystem::path pdb_path, out_dir, log_dir;
};

std::optional<arguments> read_args(int argc, const char* argv[]);

std::string app_full_name();

std::string joinStrings(std::vector<std::string> const& values, std::string const& delimiter);

