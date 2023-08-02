#pragma once

#include <iostream>
#include <fstream>
#include <vector>

class csvfile
{
    std::fstream fs;
    const std::string separator;
    bool newrow;

public:
    csvfile(const std::string& filename_, const std::string& separator_, const std::vector<std::string>& column_names);
    ~csvfile();

    void flush();
    void endrow();

    template <typename T>
    csvfile& operator<<(const T& x);
};
