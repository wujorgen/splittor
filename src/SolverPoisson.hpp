#ifndef SOLVER_POISSON
#define SOLVER_POISSON

#include <Eigen/Dense>

void SolvePoissonGrid2D(Eigen::ArrayXXd& phi, const Eigen::ArrayXXd& b, const double dx, const double dy, const double etol = 1e-4, const double MAXITERP = 100);

// void SolvePoissonGrid3D();
// Eigen does not natively support ArrayXXXd. Would have to either switch to arrayfire, or use a vector of ArrayXXd's.
// I'm fine with only running 2D stuff for now, so we stick with Eigen.
// For the unstructured stuff I want to do later, this type of data structure won't matter anyways.

#endif