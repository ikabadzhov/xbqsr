#include "Global.h"


CigarOperationsTable cigarOperationsTable;

/*
std::vector<char> nibbleToBase({static_cast<std::byte>('='), static_cast<std::byte>('A'),
                                     static_cast<std::byte>('C'), static_cast<std::byte>('M'),
                                     static_cast<std::byte>('G'), static_cast<std::byte>('R'),
                                     static_cast<std::byte>('S'), static_cast<std::byte>('V'),
                                     static_cast<std::byte>('T'), static_cast<std::byte>('W'),
                                     static_cast<std::byte>('Y'), static_cast<std::byte>('H'),
                                     static_cast<std::byte>('K'), static_cast<std::byte>('D'),
                                     static_cast<std::byte>('B'), static_cast<std::byte>('N')});*/


std::vector<char> nibbleToBase({'=', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'});

