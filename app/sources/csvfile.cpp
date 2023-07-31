#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <string>
#include "csvfile.h"


csvfile::csvfile(const std::string& filename, const std::string& separator_, const std::vector<std::string>& column_names) : fs(), separator(separator_), newrow(true)
{
    fs.exceptions(std::ios::failbit | std::ios::badbit);
    fs.open(filename);
    for (const auto& s : column_names)
        *this << s;
    endrow();
}

csvfile::~csvfile()
{
    flush();
    fs.close();
}

void csvfile::flush()
{
    fs.flush();
}

void csvfile::endrow()
{
    fs << std::endl;
    newrow = true;
}

template <class T>
csvfile& csvfile::operator<<(const T& x)
{
    std::ostringstream oss;
    oss << x;
    std::string s(oss.str());
    if (!newrow) fs << separator;
    newrow = false;
    if (s.find(separator) != std::string::npos) fs << std::quoted(s);
    else fs << s;
    return *this;
}
