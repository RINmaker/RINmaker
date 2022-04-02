#pragma once

#include <iostream>
#include <string>
#include <optional>
#include <filesystem>

#include <spdlog/spdlog.h>

#include "rin_params.h"
#include "log_manager.h"

struct arguments final
{
    rin::parameters params;
    std::filesystem::path pdb_path, out_path, log_path;
};

bool read_args(int argc, const char* argv[], std::optional<arguments>& result);

std::string app_full_name();

/*
// TODO testare che funzioni, presa da stackoverflow
std::string sprintf(const std::string& fmt, ...)
{
    int size = 100;
    std::string str;
    va_list ap;

    while (1)
    {
        str.resize(size);
        va_start(ap, fmt);
        int n = vsnprintf(&str[0], size, fmt.c_str(), ap);
        va_end(ap);

        if (n > -1 && n < size)
        {
            str.resize(n); // Make sure there are no trailing zero char
            return str;
        }
        if (n > -1)
            size = n + 1;
        else
            size *= 2;
    }
}
*/