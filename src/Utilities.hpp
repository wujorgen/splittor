#ifndef UTILITIES
#define UTILITIES

#include "Types.hpp"

#include <Eigen/Dense>

void applyVelocityBoundaryConditions2D(Eigen::VectorXd& u, Eigen::VectorXd& v, const BoundaryConditions& BC, const GridInfo& Grid);

Eigen::VectorXd ConvertFieldToVector(const Eigen::ArrayXXd& field);

Eigen::ArrayXXd ConvertVectorToField(const Eigen::VectorXd& vector, const int& NX, const int& NY);

double computeMaxDivergence(const Eigen::VectorXd& u, const Eigen::VectorXd& v, const GridInfo& Grid, const bool& print = false);

double computeMeanDivergence(const Eigen::VectorXd& u, const Eigen::VectorXd& v, const GridInfo& Grid, const bool& print = false);

#endif