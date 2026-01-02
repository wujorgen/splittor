#ifndef VELOCITY_PROJECTION
#define VELOCITY_PROJECTION

#include <Eigen/Dense>

#include "Types.hpp"

void projectVelocity(Eigen::VectorXd& u_next, Eigen::VectorXd& v_next,
    const Eigen::VectorXd& u_star, const Eigen::VectorXd& v_star, const Eigen::VectorXd& p,
    const GridInfo& Grid, const BoundaryConditions& BC, const ProblemInformation& Problem);

#endif