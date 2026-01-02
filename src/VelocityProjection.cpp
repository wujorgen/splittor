#include <Eigen/Dense>

#include "Types.hpp"
#include "VelocityProjection.hpp"

/**
 * @brief Corrects velocity field using gradient of pressure field.
 *
 * @param u_next
 * @param v_next
 * @param u_star
 * @param v_star
 * @param p
 * @param Grid
 * @param BC
 * @param Problem
 */
void projectVelocity(Eigen::VectorXd& u_next, Eigen::VectorXd& v_next,
    const Eigen::VectorXd& u_star, const Eigen::VectorXd& v_star, const Eigen::VectorXd& p,
    const GridInfo& Grid, const BoundaryConditions& BC, const ProblemInformation& Problem)
{
    // indexing
    int ij;
    int ip1, im1; // (i+1,j) and (i-1,j)
    int jp1, jm1; // (i,j+1) and (i,j-1)
    // grid steps
    double dx_p, dx_m;
    double dy_p, dy_m;
    // pressure derivatives
    double dpdx_p, dpdx_m;
    double dpdy_p, dpdy_m;
    //
    for (int jdx = 0; jdx < Grid.NY; jdx++) {
        for (int idx = 0; idx < Grid.NX; idx++) {
            // indexing
            ij = jdx * Grid.NX + idx;
            im1 = jdx * Grid.NX + idx - 1;
            ip1 = jdx * Grid.NX + idx + 1;
            jm1 = (jdx - 1) * Grid.NX + idx;
            jp1 = (jdx + 1) * Grid.NX + idx;
            //
            if (BC.velocity_type[ij] == BoundaryConditionType::DIRECHLET) {
                u_next(ij) = BC.u[ij];
                v_next(ij) = BC.v[ij];
                continue;
            }
            // set grid steps
            dx_p = Grid.X(idx + 1) - Grid.X(idx);
            dx_m = Grid.X(idx) - Grid.X(idx - 1);
            dy_p = Grid.Y(jdx + 1) - Grid.Y(jdx);
            dy_m = Grid.Y(jdx) - Grid.Y(jdx - 1);
            // calc pressure derivatives
            dpdx_p = (p(ip1) - p(ij)) / dx_p;
            dpdx_m = (p(ij) - p(im1)) / dx_m;
            dpdy_p = (p(jp1) - p(ij)) / dx_p;
            dpdy_m = (p(ij) - p(jm1)) / dy_m;
            //
            u_next(ij) = u_star(ij) - (Problem.dt / Problem.Properties.rho) * (dpdx_p + dpdx_m) / 2.0;
            v_next(ij) = v_star(ij) - (Problem.dt / Problem.Properties.rho) * (dpdy_p + dpdy_m) / 2.0;
        }
    }
}