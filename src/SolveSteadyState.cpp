#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "Constants.hpp"
#include "PISO.hpp"
#include "SolveSteadyState.hpp"
#include "Types.hpp"
#include "Utilities.hpp"

/**
 * @brief
 *
 * @param
 */
void solveSteadyStateProblem(Eigen::VectorXd& u, Eigen::VectorXd& v, Eigen::VectorXd& p,
    const GridInfo& Grid, const BoundaryConditions& BC, const ProblemInformation& Problem)
{
    // std::ofstream logfile(Constants::FNAME_LOG, std::ios::app);

    Eigen::VectorXd u_next = Eigen::VectorXd::Zero(Grid.NX * Grid.NY);
    Eigen::VectorXd v_next = Eigen::VectorXd::Zero(Grid.NX * Grid.NY);
    Eigen::VectorXd p_next = Eigen::VectorXd::Zero(Grid.NX * Grid.NY);

    double delta_u;
    double delta_v;
    std::vector<double> delta_u_history;
    std::vector<double> delta_v_history;

    int FLAG_PISO;

    for (int istep = 0; istep < Problem.Convergence.Steady.MAX_STEPS; istep++) {
        std::cout << "Begin Steady State Step " << istep << std::endl;
        applyVelocityBoundaryConditions2D(u, v, BC, Grid);
        FLAG_PISO = stepPISO(u_next, v_next, p_next, u, v, p, Grid, BC, Problem);
        if (FLAG_PISO != 0) {
            // std::cout << "PISO step failed while solving steady state problem." << std::endl;
        }
        //
        delta_u = (u_next - u).norm() / u.norm();
        delta_v = (v_next - v).norm() / v.norm();
        delta_u_history.push_back(delta_u);
        delta_v_history.push_back(delta_v);
        //
        std::cout << std::setprecision(5) << std::setfill(' ')
                  << "delta_u: " << delta_u << "\n"
                  << "delta_v: " << delta_v
                  << std::endl;
        std::cout << "End Steady State Step " << istep << "\n"
                  << std::endl;
        //
        if (delta_u < Problem.Convergence.Steady.TOL_U_REL && delta_v < Problem.Convergence.Steady.TOL_V_REL) {
            break;
        }
        //
        u = u_next;
        v = v_next;
        p = p_next;
    }
}

/**
 *
 */
void writeSteadyStateSolution(const Eigen::VectorXd& u, const Eigen::VectorXd& v, const Eigen::VectorXd& p,
    const GridInfo& Grid, const BoundaryConditions& BC, const ProblemInformation& Problem)
{
}