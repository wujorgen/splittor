#ifndef INPUTPROCESSOR
#define INPUTPROCESSOR

#include <vector>
#include "Types.hpp"

void applyBoundaryConditions2D(Eigen::VectorXd& u, Eigen::VectorXd& v, const BoundaryConditions& BC, const GridInfo& Grid);

void calcGridSteps(GridInfo &Grid);

bool checkGridIsCartesian(const GridInfo &Grid, const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z, double tolerance = 1e-6);

#endif