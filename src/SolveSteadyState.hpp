#ifndef SOLVE_STEADY_STATE
#define SOLVE_STEADY_STATE

#include <Eigen/Dense>

int solveSteadyStateProblem(Eigen::VectorXd& u_final, Eigen::VectorXd& v_final, Eigen::VectorXd& p_final,
    const GridInfo& Grid, const BoundaryConditions& BC, const ProblemInformation& Problem);

void writeSteadyStateSolution(Eigen::VectorXd& u, Eigen::VectorXd& v, Eigen::VectorXd& p,
    const GridInfo& Grid, const BoundaryConditions& BC, const ProblemInformation& Problem);

#endif