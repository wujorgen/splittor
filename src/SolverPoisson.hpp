#ifndef SOLVER_POISSON
#define SOLVER_POISSON

#include <Eigen/Dense>

void SolvePoisson(Eigen::ArrayXXd& phi, const Eigen::ArrayXXd& b, double etol = 1e-4);

#endif