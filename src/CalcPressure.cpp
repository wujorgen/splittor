#include <Eigen/Dense>

#include "CalcPressure.hpp"
#include "Types.hpp"

/**
 * @brief 
 */
void calcPressure(Eigen::VectorXd& p_star,
    const Eigen::VectorXd& u_star, const Eigen::VectorXd& v_star, 
    const GridInfo& Grid, const BoundaryConditions& BC, const ProblemInformation& Problem)
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
    double dvdy_p, dvdy_m;
    // coefficient matrix and RHS
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(Grid.NX * Grid.NY, Grid.NX * Grid.NY);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(Grid.NX * Grid.NY);

    for (int jdx = 0; jdx < Grid.NY; jdx++) {
        for (int idx = 0; idx < Grid.NX; idx++) {
            ij = jdx * Grid.NX + idx;
            im1 = jdx * Grid.NX + idx - 1;
            ip1 = jdx * Grid.NX + idx + 1;
            jm1 = (jdx - 1) * Grid.NX + idx;
            jp1 = (jdx + 1) * Grid.NX + idx;
            // handle boundaries
            // by default, we assume no pressure at the boundaries of the domain unless otherwise specified.
            if (BC.pressure_type[ij] == BoundaryConditionType::DIRECHLET) {
                A(ij, ij) = 1;
                b(ij) = 0;
                continue;
            } else if (jdx == Grid.NY - 1) {
                A(ij, ij) = 1;
                A(ij, jm1) = -1;
                b(ij) = 0;
                continue;
            } else if (jdx == 0) {
                A(ij, ij) = 1;
                A(ij, jp1) = -1;
                b(ij) = 0;
                continue;
            } else if (idx == Grid.NX - 1) {
                A(ij, ij) = 1;
                A(ij, im1) = -1;
                b(ij) = 0;
                continue;
            } else if (idx == 0) {
                A(ij, ij) = 1;
                A(ij, ip1) = -1;
                b(ij) = 0;
                continue;
            }
            // set grid steps
            dx_p = Grid.X(idx + 1) - Grid.X(idx);
            dx_m = Grid.X(idx) - Grid.X(idx - 1);
            dy_p = Grid.Y(jdx + 1) - Grid.Y(jdx);
            dy_m = Grid.Y(jdx) - Grid.Y(jdx - 1);
            // calculate first order derivatives
            dudx_p = (u_star(ip1) - u_star(ij)) / dx_p;
            dudx_m = (u_star(ij) - u_star(im1)) / dx_m;
            dvdy_p = (v_star(jp1) - v_star(ij)) / dy_p;
            dvdy_m = (v_star(ij) - v_star(jm1)) / dy_m;
            // build coefficient matrix
            A(ij, ij) = -2.0 * (dx_p + dx_m) / (dx_p * dx_m * (dx_p + dx_m)) - 2.0 * (dy_p + dy_m) / (dy_p * dy_m * (dy_p + dy_m));
            A(ij, ip1) = 2.0 * dx_m / (dx_p * dx_m * (dx_p + dx_m));
            A(ij, im1) = 2.0 * dx_p / (dx_p * dx_m * (dx_p + dx_m));
            A(ij, jp1) = 2.0 * dy_m / (dy_p * dy_m * (dy_p + dy_m));
            A(ij, jm1) = 2.0 * dy_p / (dy_p * dy_m * (dy_p + dy_m));
            // build RHS
            b(ij) = Problem.Properties.rho / Problem.dt * (0.5 * (dudx_p + dudx_m) + 0.5 * (dvdy_p + dvdy_m));
        }
    }
    // solve
    p_star = A.partialPivLu().solve(b);
}