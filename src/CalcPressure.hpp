#ifndef CALC_PRESSURE
#define CALC_PRESSURE

#include <Eigen/Dense>

#include "Types.hpp"

void calcPressure(Eigen::VectorXd& p_star,
    const Eigen::VectorXd& u_star, const Eigen::VectorXd& v_star,
    const GridInfo& Grid, const BoundaryConditions& BC, const ProblemInformation& Problem);

#endif