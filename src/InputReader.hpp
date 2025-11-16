#ifndef INPUTREADER
#define INPUTREADER

#include "Types.hpp"

int ReadInputFile(ProblemInformation&, GridInfo&, BoundaryConditions&, const std::string& fname = "grid.txt");

#endif