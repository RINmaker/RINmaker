#pragma once

#include <iostream>

#include <string>
#include <vector>

#include <algorithm>
#include <functional>

#include <cctype>
#include <locale>

#include <set>

#include <iostream>
#include <string>
#include <memory>
#include <cstdio>

#include "CLI/CLI.hpp"

#include "config.h"

// #include <varargs.h>

namespace prelude {
// trim from start (in place)
// https://stackoverflow.com/a/217605
    inline static void ltrim(std::string &s) {
        s.erase(s.begin(), find_if(s.begin(), s.end(), [](unsigned char ch) { return !isspace(ch); }));
    }

// trim from both ends
// https://stackoverflow.com/a/217605
    inline static void rtrim(std::string &s) {
        s.erase(find_if(s.rbegin(), s.rend(), [](unsigned char ch) { return !isspace(ch); }).base(), s.end());
    }

    inline static std::string trim(std::string s) {
        ltrim(s);
        rtrim(s);
        return s;
    }

    inline static std::string sort(std::string const &id_1, std::string const &id_2) {
        return id_1 < id_2 ? id_1 + id_2 : id_2 + id_1;
    }

    inline static bool match(std::string const &str, std::vector<std::string> const &patterns) {
        for (auto const &p: patterns)
            if (str.find(p) != std::string::npos)
                return true;

        return false;
    }

    inline static bool match(std::string const &str, std::string const &pattern) {
        return str.find(pattern) != std::string::npos;
    }
}
