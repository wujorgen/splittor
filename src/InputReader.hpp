#ifndef INPUTREADER
#define INPUTREADER

#include "Types.hpp"

int readGridFile(GridInfo& Grid, BoundaryConditions& BC, const std::string& fname = "grid.txt");

int readProblemInformation(ProblemInformation& Problem, const std::string& fname = "probleminformation.txt");

#endif