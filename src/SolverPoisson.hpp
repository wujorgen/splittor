#ifndef SOLVER_POISSON
#define SOLVER_POISSON

#include <Eigen/Dense>

void SolvePoisson(Eigen::ArrayXXd& phi, const Eigen::ArrayXXd& b, const double dx, const double dy, const double etol = 1e-4, const double MAXITERP = 100);

#endif