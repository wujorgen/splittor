#ifndef INPUTREADER
#define INPUTREADER

#include "Types.hpp"

int ReadGridFile(ProblemInformation& Problem, GridInfo& Grid, BoundaryConditions& BC, const std::string& fname = "grid.txt");

int ReadProblemInformation(ProblemInformation& Problem, GridInfo& Grid, BoundaryConditions& BC, const std::string& fname = "probleminfo.txt");

#endif