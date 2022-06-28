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
}
