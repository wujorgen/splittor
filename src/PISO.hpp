#ifndef PISO
#define PISO

#include <Eigen/Dense>
#include <Types.hpp>

void stepPISO(Eigen::VectorXd& u_next, Eigen::VectorXd& v_next, Eigen::VectorXd& p_next,
    const Eigen::VectorXd& u, const Eigen::VectorXd& v, const Eigen::VectorXd& p,
    const GridInfo& Grid, const BoundaryConditions& BC, const ProblemInformation& Problem);

#endif