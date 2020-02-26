#include "parse.h"

#include <fstream>
#include <iostream>


void parseCsv()
{
    std::ifstream ifs("in/est_states.floats", std::ios::binary);
    std::ifstream::pos_type pos = ifs.tellg();

    float f;
    while (ifs.read(reinterpret_cast<char*>(&f), sizeof(float))) std::cout << f << '\n';
}

