#ifndef UTILITIES
#define UTILITIES

#include "Types.hpp"

#include <Eigen/Dense>

void applyVelocityBoundaryConditions2D(Eigen::VectorXd& u, Eigen::VectorXd& v, const BoundaryConditions& BC, const GridInfo& Grid);

Eigen::VectorXd ConvertFieldToVector(const Eigen::ArrayXXd& field);

Eigen::ArrayXXd ConvertVectorToField(const Eigen::VectorXd& vector, const int& NX, const int& NY);

#endif