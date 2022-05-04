#pragma once

#include <string>
#include <array>

#include "prelude.h"

namespace record
{
template<size_t N, typename Derived>
class base
{
protected:
    uint32_t _line_number;
    std::array<std::string, N> _fields;

    explicit base(std::string const& line, uint32_t line_number) : _line_number(line_number)
    {
        size_t line_size = line.size();
        for (size_t i = 0; i < N; ++i)
        {
            size_t first = Derived::_limits[i].first;
            size_t second = Derived::_limits[i].second;

            if (first + second <= line_size)
                _fields[i] = prelude::trim(line.substr(first, second));
        }
    }

public:
    [[nodiscard]]
    uint32_t line_number() const
    { return _line_number; }
};

class atom : public base<15, atom>
{
private:
    friend class base<15, atom>;

    // https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
    constexpr static std::array<std::pair<size_t, size_t>, 15> const _limits = {
            std::make_pair(0, 6),   //  0 "ATOM"
            std::make_pair(6, 5),   //  1 serial
            std::make_pair(12, 4),  //  2 name
            std::make_pair(16, 1),  //  3 altLoc
            std::make_pair(17, 3),  //  4 resName
            std::make_pair(21, 1),  //  5 chainID
            std::make_pair(22, 4),  //  6 resSeq
            std::make_pair(26, 1),  //  7 iCode
            std::make_pair(30, 8),  //  8 x
            std::make_pair(38, 8),  //  9 y
            std::make_pair(46, 8),  // 10 z
            std::make_pair(54, 6),  // 11 occupancy
            std::make_pair(60, 6),  // 12 tempFactor
            std::make_pair(76, 2),  // 13 element
            std::make_pair(78, 2)   // 14 charge
    };

public:
    explicit atom(std::string const& line, uint32_t line_number) : base(line, line_number)
    {}

    [[nodiscard]]
    bool same_res(atom const& other) const
    {
        return chain_id() == other.chain_id() && res_name() == other.res_name() && res_seq() == other.res_seq();
    }

    [[nodiscard]]
    std::string const& name() const
    { return _fields[2]; }

    [[nodiscard]]
    std::string const& res_name() const
    { return _fields[4]; }

    [[nodiscard]]
    std::string const& chain_id() const
    { return _fields[5]; }

    [[nodiscard]]
    std::string const& element_name() const
    { return _fields[13]; }

    [[nodiscard]]
    int serial() const
    { return std::stoi(_fields[1]); }

    [[nodiscard]]
    int charge() const
    {
        auto const& x = _fields[14];
        return x == "1+" ? 1 : x == "1-" ? -1 : 0;
    } // @TODO config

    [[nodiscard]]
    int res_seq() const
    { return std::stoi(_fields[6]); }

    [[nodiscard]]
    double x() const
    { return std::stod(_fields[8]); }

    [[nodiscard]]
    double y() const
    { return std::stod(_fields[9]); }

    [[nodiscard]]
    double z() const
    { return std::stod(_fields[10]); }

    [[nodiscard]]
    double temp_factor() const
    { return std::stod(_fields[12]); }
};

class helix : public base<14, helix>
{
private:
    friend class base<14, helix>;

    // http://www.wwpdb.org/documentation/file-format-content/format33/sect5.html#HELIX
    constexpr static std::array<std::pair<size_t, size_t>, 14> const _limits = {
            std::make_pair<size_t, size_t>(0, 6),   // HELIX
            std::make_pair<size_t, size_t>(7, 3),   // serial number (starts from 1)
            std::make_pair<size_t, size_t>(11, 3),  // helix ID
            std::make_pair<size_t, size_t>(15, 3),  // res name of initial residue
            std::make_pair<size_t, size_t>(19, 1),  // chain id of initial residue
            std::make_pair<size_t, size_t>(21, 4),  // seq number of initial residue
            std::make_pair<size_t, size_t>(25, 1),  // insertion code of initial residue
            std::make_pair<size_t, size_t>(27, 3),  // res name of final residue
            std::make_pair<size_t, size_t>(31, 1),  // chain id of final residue
            std::make_pair<size_t, size_t>(33, 4),  // eq number of final residue
            std::make_pair<size_t, size_t>(37, 1),  // insertion code of final residue
            std::make_pair<size_t, size_t>(38, 1),  // helix class
            std::make_pair<size_t, size_t>(40, 30), // helix comment
            std::make_pair<size_t, size_t>(71, 5)   // helix length
    };

public:
    explicit helix(std::string const& line, uint32_t line_number) : base(line, line_number)
    {}

    [[nodiscard]]
    std::string const& init_chain_id() const
    { return _fields[4]; }

    [[nodiscard]]
    std::string const& init_res_name() const
    { return _fields[3]; }

    [[nodiscard]]
    int init_seq_number() const
    { return std::stoi(_fields[5]); }

    [[nodiscard]]
    std::string const& end_chain_id() const
    { return _fields[8]; }

    [[nodiscard]]
    std::string const& end_res_name() const
    { return _fields[7]; }

    [[nodiscard]]
    int end_seq_number() const
    { return std::stoi(_fields[9]); }

    [[nodiscard]]
    std::string get_id() const
    { return _fields[2]; }

    [[nodiscard]]
    std::string init_res_id() const
    {
        return init_chain_id() + ":" + std::to_string(init_seq_number()) + ":_:" + init_res_name();
    }

    [[nodiscard]]
    std::string end_res_id() const
    {
        return end_chain_id() + ":" + std::to_string(end_seq_number()) + ":_:" + end_res_name();
    }

    [[nodiscard]]
    prelude::interval<int> range() const
    { return {init_seq_number(), end_seq_number()}; }

    [[nodiscard]]
    int serial() const
    { return std::stoi(_fields[1]); }
};

class sheet_piece : public base<23, sheet_piece>
{
private:
    friend class base<23, sheet_piece>;

    // http://www.wwpdb.org/documentation/file-format-content/format33/sect5.html#SHEET
    constexpr static std::array<std::pair<size_t, size_t>, 23> const _limits = {
            std::make_pair<size_t, size_t>(0, 6),    // SHEET
            std::make_pair<size_t, size_t>(7, 3),    // strand number (starts with 0)
            std::make_pair<size_t, size_t>(11, 3),   // sheet id
            std::make_pair<size_t, size_t>(14, 2),   // number of strand in the sheet
            std::make_pair<size_t, size_t>(17, 3),   // name of initial residue
            std::make_pair<size_t, size_t>(21, 1),   // chain id of initial residue
            std::make_pair<size_t, size_t>(22, 4),   // seq number of initial residue
            std::make_pair<size_t, size_t>(26, 1),   // insertion code of initial residue
            std::make_pair<size_t, size_t>(28, 3),   // nome of final residue
            std::make_pair<size_t, size_t>(32, 1),   // chain id of final residue
            std::make_pair<size_t, size_t>(33, 4),   // seq number of final residue
            std::make_pair<size_t, size_t>(37, 1),   // insertion code of final residue
            std::make_pair<size_t, size_t>(38, 2),   // sense of strand
            std::make_pair<size_t, size_t>(41, 4),   // atom name of current strand
            std::make_pair<size_t, size_t>(45, 3),   // atom res name of current strand
            std::make_pair<size_t, size_t>(49, 1),   // atom chain id of current strand
            std::make_pair<size_t, size_t>(50, 4),   // atom res seq of current strand
            std::make_pair<size_t, size_t>(54, 1),   // atom insertion code of current strand
            std::make_pair<size_t, size_t>(56, 4),   // atom name of precedent strand
            std::make_pair<size_t, size_t>(60, 3),   // atom res name of precedent strand
            std::make_pair<size_t, size_t>(64, 1),   // atom chain id of precedent strand
            std::make_pair<size_t, size_t>(65, 4),   // atom res seq of precedent strand
            std::make_pair<size_t, size_t>(69, 1)    // atom insertion code of precedent strand
    };

public:
    explicit sheet_piece(std::string const& line, uint32_t line_number) : base(line, line_number)
    {}

    [[nodiscard]]
    std::string const& init_chain_id() const
    { return _fields[5]; }

    [[nodiscard]]
    std::string const& init_res_name() const
    { return _fields[4]; }

    [[nodiscard]]
    int init_seq_number() const
    { return std::stoi(_fields[6]); }

    [[nodiscard]]
    std::string const& end_chain_id() const
    { return _fields[9]; }

    [[nodiscard]]
    std::string const& end_res_name() const
    { return _fields[8]; }

    [[nodiscard]]
    int end_seq_number() const
    { return std::stoi(_fields[10]); }

    [[nodiscard]]
    std::string get_id() const
    { return _fields[2]; }

    [[nodiscard]]
    std::string init_res_id() const
    {
        return init_chain_id() + ":" + std::to_string(init_seq_number()) + ":_:" + init_res_name();
    }

    [[nodiscard]]
    std::string end_res_id() const
    {
        return end_chain_id() + ":" + std::to_string(end_seq_number()) + ":_:" + end_res_name();
    }

    [[nodiscard]]
    prelude::interval<int> range() const
    { return {init_seq_number(), end_seq_number()}; }

    [[nodiscard]]
    int incremental_strand_number() const
    { return std::stoi(_fields[1]); }
};

class ss : public base<13, ss>
{
private:
    friend class base<13, ss>;

    // https://www.wwpdb.org/documentation/file-format-content/format33/sect6.html#SSBOND
    constexpr static std::array<std::pair<size_t, size_t>, 13> const _limits = {
            std::make_pair(0, 6),   // "SSBOND"
            std::make_pair(7, 3),   // serNum
            std::make_pair(11, 3),  // "CYS"
            std::make_pair(15, 1),  // chainID1
            std::make_pair(17, 4),  // seqNum1
            std::make_pair(21, 1),  // icode1
            std::make_pair(25, 3),  // "CYS"
            std::make_pair(29, 1),  // chainID2
            std::make_pair(31, 4),  // seqNum2
            std::make_pair(35, 1),  // icode2
            std::make_pair(59, 6),  // sym1
            std::make_pair(66, 6),  // sym2
            std::make_pair(73, 5)   // length
    };

public:
    explicit ss(std::string const& line, uint32_t line_number) : base(line, line_number)
    {}

    [[nodiscard]]
    double length() const
    { return std::stod(_fields[12]); }

    [[nodiscard]]
    int seq_num_1() const
    { return std::stoi(_fields[4]); }

    [[nodiscard]]
    int seq_num_2() const
    { return std::stoi(_fields[8]); }

    [[nodiscard]]
    std::string const& name_1() const
    { return _fields[2]; }

    [[nodiscard]]
    std::string const& name_2() const
    { return _fields[6]; }

    [[nodiscard]]
    std::string const& chain_id_1() const
    { return _fields[3]; }

    [[nodiscard]]
    std::string const& chain_id_2() const
    { return _fields[7]; }
};
}
