#ifndef INPUTPROCESSOR
#define INPUTPROCESSOR

#include <vector>
#include "Types.hpp"

void CalcGridSteps(GridInfo& Grid);

bool GridIsCartesian(const GridInfo& Grid, const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, double tolerance = 1e-6);

#endif