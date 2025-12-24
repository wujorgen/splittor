#include <Eigen/Dense>

#include "CalcIntermediateVelocity.hpp"
#include "Types.hpp"

/**
 * @brief Calculates intermediate velocities using an explicit timestep.
 *
 * As this function utilizes an explcit step, it is subject to CFL number for stability.
 * @param u_star
 * @param v_star
 * @param u
 * @param v
 * @param Grid
 * @param BC
 */
void stepIntermediateExplicit(Eigen::VectorXd& u_star, Eigen::VectorXd& v_star, const Eigen::VectorXd& u, const Eigen::VectorXd& v, const GridInfo& Grid, const BoundaryConditions& BC, const ProblemInformation& Problem)
{
    // indexing
    int ij;
    int ip1, im1; // (i+1,j) and (i-1,j)
    int jp1, jm1; // (i,j+1) and (i,j-1)
    // grid steps
    double dx_p, dx_m;
    double dy_p, dy_m;
    // first order derivatives
    double dudx_p, dudx_m;
    double dudy_p, dudy_m;
    double dvdx_p, dvdx_m;
    double dvdy_p, dvdy_m;
    // nonlinear first order terms
    double ududx, vdudy;
    double udvdx, vdvdy;
    // second order terms
    double d2udx2, d2udy2;
    double d2vdx2, d2vdy2;

    for (int jdx = 0; jdx < Grid.NY; jdx++) {
        for (int idx = 0; idx < Grid.NX; idx++) {
            ij = jdx * Grid.NX + idx;
            im1 = jdx * Grid.NX + idx - 1;
            ip1 = jdx * Grid.NX + idx + 1;
            jm1 = (jdx - 1) * Grid.NX + idx;
            jp1 = (jdx + 1) * Grid.NX + idx;
            // handle boundary conditions
            if (false) {
                // TODO: depends on contents of BC
                continue;
            }
            // set grid steps
            dx_p = Grid.dx(idx + 1);
            dx_m = Grid.dx(idx);
            dy_p = Grid.dy(jdx + 1);
            dy_m = Grid.dy(jdx);
            // calculate first order derivatives
            dudx_p = (u(ip1) - u(ij)) / dx_p;
            dudx_m = (u(ij) - u(im1)) / dx_m;
            dudy_p = (u(jp1) - u(ij)) / dy_p;
            dudy_m = (u(ij) - u(jm1)) / dy_m;
            dvdx_p = (v(ip1) - v(ij)) / dx_p;
            dvdx_m = (v(ij) - v(im1)) / dx_m;
            dvdy_p = (v(jp1) - v(ij)) / dy_p;
            dvdy_m = (v(ij) - v(jm1)) / dy_m;
            // calculate nonlinear first order terms
            ududx = u(ij) * (dudx_p + dudx_m) / 2;
            vdudy = v(ij) * (dudy_p + dudy_m) / 2;
            udvdx = u(ij) * (dvdx_p + dvdx_m) / 2;
            vdvdy = v(ij) * (dvdy_p + dvdy_m) / 2;
            // calculate second order terms
            d2udx2 = (dudx_p - dudx_m) / ((dx_p + dx_m) / 2);
            d2udy2 = (dudy_p - dudy_m) / ((dy_p + dy_m) / 2);
            d2vdx2 = (dvdx_p - dvdx_m) / ((dx_p + dx_m) / 2);
            d2vdy2 = (dvdy_p - dvdy_m) / ((dy_p + dy_m) / 2);
            // calculate intermediate velocity
            u_star(ij) = u(ij) + Problem.dt * (-ududx - vdudy + (Problem.Properties.mu / Problem.Properties.rho) * (d2udx2 + d2udy2));
            v_star(ij) = v(ij) + Problem.dt * (-udvdx - vdvdy + (Problem.Properties.mu / Problem.Properties.rho) * (d2vdx2 + d2vdy2));
        }
    }
}

/**
 * @brief Calculates intermediate velocities using a semi implicit timestep.
 * 
 * The LHS of the equation contains the timestep and diffusion terms, while the nonlinear advetion term is moved to the RHS.
 * This eliminates the need for a nonlinear implicit solve.
 * @param u_star
 * @param v_star
 * @param u
 * @param v
 * @param Grid
 * @param BC
 */
void stepIntermediateSemiImplicit(Eigen::VectorXd& u_star, Eigen::VectorXd& v_star, const Eigen::VectorXd& u, const Eigen::VectorXd& v, const GridInfo& Grid, const BoundaryConditions& BC, const ProblemInformation& Problem)
{
    // indexing
    int ij;
    int ip1, im1; // (i+1,j) and (i-1,j)
    int jp1, jm1; // (i,j+1) and (i,j-1)
}