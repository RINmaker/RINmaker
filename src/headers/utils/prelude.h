#pragma once

#include <iostream>

#include <string>
#include <vector>

#include <algorithm>
#include <functional>

#include <cctype>
#include <locale>

#include <set>

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
