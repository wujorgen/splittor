#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "CalcIntermediateVelocity.hpp"
#include "CalcPressure.hpp"
#include "Constants.hpp"
#include "Types.hpp"
#include "Utilities.hpp"
#include "VelocityProjection.hpp"

/**
 * @brief Calculates the next velocity and pressure fields according ot the PISO algorithm
 *
 * @param u_next
 * @param v_next
 * @param p_next
 * @param u
 * @param v
 * @param p
 * @param Grid
 * @param BC
 * @param Problem
 * @return 0 if converged solution has been stored to u_next, v_next, p_next
 */
int stepPISO(Eigen::VectorXd& u_next, Eigen::VectorXd& v_next, Eigen::VectorXd& p_next,
    const Eigen::VectorXd& u, const Eigen::VectorXd& v, const Eigen::VectorXd& p,
    const GridInfo& Grid, const BoundaryConditions& BC, const ProblemInformation& Problem)
{
    // std::ofstream logfile(Constants::FNAME_LOG, std::ios::app);

    int RETURN_FLAG = 1;
    bool DIVERGENCE_CONVERGED = false;
    bool VELOCITY_CONVERGED = false;

    double step_divergence;
    std::vector<double> divergence_history;
    double delta_u_rel_norm;
    double delta_v_rel_norm;
    std::vector<double> delta_u_history;
    std::vector<double> delta_v_history;

    Eigen::VectorXd u_s1 = Eigen::VectorXd::Zero(Grid.NX * Grid.NY);
    Eigen::VectorXd v_s1 = Eigen::VectorXd::Zero(Grid.NX * Grid.NY);
    Eigen::VectorXd p_s1 = Eigen::VectorXd::Zero(Grid.NX * Grid.NY);
    Eigen::VectorXd u_s2 = Eigen::VectorXd::Zero(Grid.NX * Grid.NY);
    Eigen::VectorXd v_s2 = Eigen::VectorXd::Zero(Grid.NX * Grid.NY);

    // solve momentum for intermediate velocity
    stepIntermediateSemiImplicit(u_s1, v_s1, u, v, Grid, BC, Problem);

    // a projection step results in a divergence free velocity field
    // however, the most last calculated velocity field and last calculated pressure field may not satisfy the pressure poisson equation
    // as such, iteration is required.
    std::cout << "Begin PISO Step" << std::endl;
    std::cout << std::setfill(' ') << std::setw(2) << "iteration | divergence | u rel norm | v rel norm" << std::endl;
    for (int itr = 0; itr < Problem.Convergence.PISO.MAX_NUM_ITERATIONS; itr++) {
        // Divergence
        step_divergence = computeDivergence(u_s1, v_s1, Grid);
        if (step_divergence < Problem.Convergence.PISO.TOL_DIVERGENCE) {
            DIVERGENCE_CONVERGED = true;
        }
        divergence_history.push_back(step_divergence);
        // step PISO
        calcPressure(p_s1, u_s1, v_s1, Grid, BC, Problem);
        projectVelocity(u_s2, v_s2, u_s1, v_s1, p_s1, Grid, BC, Problem);
        // velocity changes
        delta_u_rel_norm = (u_s2 - u_s1).norm() / u_s1.norm();
        delta_v_rel_norm = (v_s2 - v_s1).norm() / v_s1.norm();
        if (delta_u_rel_norm < Problem.Convergence.PISO.TOL_U_REL && delta_v_rel_norm < Problem.Convergence.PISO.TOL_V_REL) {
            VELOCITY_CONVERGED = true;
        }
        delta_u_history.push_back(delta_u_rel_norm);
        delta_v_history.push_back(delta_v_rel_norm);
        // print formating
        if (itr % Problem.Convergence.PISO.PRINT_EVERY_N_ITERATIONS == 0 || (DIVERGENCE_CONVERGED || VELOCITY_CONVERGED)) {
            std::cout << std::setprecision(5) << std::setfill(' ')
                      << std::setw(9) << itr << " | "
                      << std::setw(10) << step_divergence << " | "
                      << std::setw(10) << delta_u_rel_norm << " | "
                      << std::setw(10) << delta_v_rel_norm
                      << std::endl;
        }
        // exit if converged
        if (DIVERGENCE_CONVERGED || VELOCITY_CONVERGED) {
            RETURN_FLAG = 0;
            break;
        }
        // prepare for next iteration
        u_s1 = u_s2;
        v_s1 = v_s2;
    }
    std::cout << "End PISO Step" << std::endl;
    // final values
    u_next = u_s2;
    v_next = v_s2;
    p_next = p_s1;
    //
    return RETURN_FLAG;
}