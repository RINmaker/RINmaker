#pragma once

#pragma warning(push, 0)

#include <CLI/CLI.hpp>

#pragma warning(pop)

#include <string>
#include <optional>
#include <variant>
#include <filesystem>

#include "rin_params.h"
#include "log_manager.h"

std::optional<rin::parameters> read_args(int argc, const char* argv[]);

std::string app_full_name();
