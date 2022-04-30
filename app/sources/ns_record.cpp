#include "ns_record.h"

using namespace records;

// FIXME: these _limits initializations should be avoided

// https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
std::array<std::pair<size_t, size_t>, 15> const atom::_limits = {
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

// http://www.wwpdb.org/documentation/file-format-content/format33/sect5.html#HELIX
std::array<std::pair<size_t, size_t>, 14> const helix::_limits = {
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

// http://www.wwpdb.org/documentation/file-format-content/format33/sect5.html#SHEET
std::array<std::pair<size_t, size_t>, 23> const sheet_piece::_limits = {
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

// https://www.wwpdb.org/documentation/file-format-content/format33/sect6.html#SSBOND
std::array<std::pair<size_t, size_t>, 13> const ss::_limits = {
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