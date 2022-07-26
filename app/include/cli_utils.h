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

struct output_file
{ std::filesystem::path value; };

struct output_directory
{ std::filesystem::path value; };

struct arguments final
{
    rin::parameters params;

    // always a file
    std::filesystem::path input;

    // might be a single file or a directory
    std::variant<output_file, output_directory> output;

    // always a directory
    std::filesystem::path log_dir;

    // speeds up things by not keeping waters
    bool skip_water;

    // execute reduce
    bool reduce;
};

std::optional<arguments> read_args(int argc, const char* argv[]);

std::string app_full_name();
