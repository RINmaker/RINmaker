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
    std::filesystem::path pdb_path, out_path, log_path;
};

std::optional<arguments> read_args(int argc, const char* argv[]);

std::vector<std::pair<uint32_t, std::string>> read_lines(std::filesystem::path const& file_path);

std::string app_full_name();

string joinStrings(std::vector<std::string> const& values, string const& delimiter);

