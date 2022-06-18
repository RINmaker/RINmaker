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

#include "config.h"

// #include <varargs.h>

namespace prelude
{

// for use in sorted structures, such as std::map
// (see std::map requirements @ cppreference.com)
//
// semantics:
// a < b means "interval a comes before interval b and they do not overlap"
// a > b means b < a
// a == b means !(a < b) && !(b < a)
//
// equality means "overlapping"
//
template<typename T>
struct interval final
{
private:
    T const _inf, _sup;

public:
    interval(T const &a, T const &b) : _inf(a < b ? a : b), _sup(a < b ? b : a)
    {}

    struct less final
    {
        bool operator()(interval const &lhs, interval const &rhs) const
        {
            return lhs._sup < rhs._inf;
        }
    };
};

// trim from start (in place)
// https://stackoverflow.com/a/217605
inline static void ltrim(std::string &s)
{
    s.erase(s.begin(), find_if(s.begin(), s.end(), [](unsigned char ch)
    { return !isspace(ch); }));
}

// trim from both ends
// https://stackoverflow.com/a/217605
inline static void rtrim(std::string &s)
{
    s.erase(find_if(s.rbegin(), s.rend(), [](unsigned char ch)
    { return !isspace(ch); }).base(), s.end());
}

inline static std::string trim(std::string s)
{
    ltrim(s);
    rtrim(s);
    return s;
}

inline static bool match(std::string const &str, std::string const &pattern)
{
    return str.find(pattern) != std::string::npos;
}
}
