#ifndef INPUTREADER
#define INPUTREADER

#include "Types.hpp"

int ReadGridFile(ProblemInformation&, GridInfo&, BoundaryConditions&, const std::string& fname = "grid.txt");

#endif