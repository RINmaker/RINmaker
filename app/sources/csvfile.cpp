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
