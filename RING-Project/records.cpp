#include "records.h"

using namespace records;

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
        std::make_pair<size_t, size_t>(7, 3),   // serial number comincia a 1 e si incrementa
        std::make_pair<size_t, size_t>(11, 3),  // helix ID
        std::make_pair<size_t, size_t>(15, 3),  // res name del residuo iniziale
        std::make_pair<size_t, size_t>(19, 1),  // chain id del resisuo inizialie
        std::make_pair<size_t, size_t>(21, 4),  // seq number del residuo iniziale
        std::make_pair<size_t, size_t>(25, 1),  // insertion code residuo iniziale
        std::make_pair<size_t, size_t>(27, 3),  // res name del residuo finale
        std::make_pair<size_t, size_t>(31, 1),  // chain id del residuo finale
        std::make_pair<size_t, size_t>(33, 4),  // eq number del residuo finale
        std::make_pair<size_t, size_t>(37, 1),  // insertion code residuo finale
        std::make_pair<size_t, size_t>(38, 1),  // helix class
        std::make_pair<size_t, size_t>(40, 30), // helix comment
        std::make_pair<size_t, size_t>(71, 5)   // lunghezza dell'elica
};

// http://www.wwpdb.org/documentation/file-format-content/format33/sect5.html#SHEET
std::array<std::pair<size_t, size_t>, 23> const sheet_piece::_limits = {
        std::make_pair<size_t, size_t>(0, 6),    // SHEET
        std::make_pair<size_t, size_t>(7, 3),    // strand number, comincia da 1 e si incrementa di 1
        std::make_pair<size_t, size_t>(11, 3),   // sheed id
        std::make_pair<size_t, size_t>(14, 2),   // numero degli strand nel sheet
        std::make_pair<size_t, size_t>(17, 3),   // name del residuo iniziale
        std::make_pair<size_t, size_t>(21, 1),   // chain id del residuo iniziale
        std::make_pair<size_t, size_t>(22, 4),   // seq number del residuo iniziale
        std::make_pair<size_t, size_t>(26, 1),   // insertion code del residuo iniziale
        std::make_pair<size_t, size_t>(28, 3),   // nome del residuo finale
        std::make_pair<size_t, size_t>(32, 1),   // chain id  del residuo finale
        std::make_pair<size_t, size_t>(33, 4),   // seq number  del residuo finale
        std::make_pair<size_t, size_t>(37, 1),   // insertion code del residuo finale
        std::make_pair<size_t, size_t>(38, 2),   // sense of strand
        std::make_pair<size_t, size_t>(41, 4),   // atom name nello strand corrente
        std::make_pair<size_t, size_t>(45, 3),   // res name dell' atomo nello strand corrente
        std::make_pair<size_t, size_t>(49, 1),   // chain id dell' atomo nello strand corrente
        std::make_pair<size_t, size_t>(50, 4),   // res seq dell' atomo nello strand corrente
        std::make_pair<size_t, size_t>(54, 1),   // insertion code dell' atomo nello strand corrente
        std::make_pair<size_t, size_t>(56, 4),   // atom name nello strand precedente
        std::make_pair<size_t, size_t>(60, 3),   // res name dell' atomo nello strand precedente
        std::make_pair<size_t, size_t>(64, 1),   // chain id dell' atomo nello strand precedente
        std::make_pair<size_t, size_t>(65, 4),   // res seq dell' atomo nello strand precedente
        std::make_pair<size_t, size_t>(69, 1)    // insertion code dell' atomo nello strand precedente
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