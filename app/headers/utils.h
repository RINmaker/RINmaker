#pragma once

#include <iostream>
#include <string>
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

bool read_args(int argc, const char* argv[], std::optional<arguments>& result);

std::string app_full_name();

string joinStrings(std::vector<std::string> const& values, string const& delimiter);

