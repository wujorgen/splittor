#ifndef CALC_INTERMEDIATE_VELOCITY
#define CALC_INTERMEDIATE_VELOCITY

#include "Types.hpp"
#include <Eigen/Dense>

void stepIntermediateExplicit(Eigen::VectorXd& u_star, Eigen::VectorXd& v_star, const Eigen::VectorXd& u, const Eigen::VectorXd& v, const GridInfo& Grid, const BoundaryConditions& BC, const ProblemInformation& Problem);

void stepIntermediateSemiImplicit(Eigen::VectorXd& u_star, Eigen::VectorXd& v_star, const Eigen::VectorXd& u, const Eigen::VectorXd& v, const GridInfo& Grid, const BoundaryConditions& BC, const ProblemInformation& Problem);

#endif