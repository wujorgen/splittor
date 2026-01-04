#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <iostream>

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
            // std::cout << "idx,jdx: " << idx << ", " << jdx << std::endl;
            ij = jdx * Grid.NX + idx;
            im1 = jdx * Grid.NX + idx - 1;
            ip1 = jdx * Grid.NX + idx + 1;
            jm1 = (jdx - 1) * Grid.NX + idx;
            jp1 = (jdx + 1) * Grid.NX + idx;
            //
            if (BC.velocity_type[ij] == BoundaryConditionType::DIRECHLET) {
                u_star(ij) = BC.u[ij];
                v_star(ij) = BC.v[ij];
                continue;
            }
            // set grid steps
            dx_p = Grid.X(idx + 1) - Grid.X(idx);
            dx_m = Grid.X(idx) - Grid.X(idx - 1);
            dy_p = Grid.Y(jdx + 1) - Grid.Y(jdx);
            dy_m = Grid.Y(jdx) - Grid.Y(jdx - 1);
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
            ududx = u(ij) * (dudx_p + dudx_m) / 2.0;
            vdudy = v(ij) * (dudy_p + dudy_m) / 2.0;
            udvdx = u(ij) * (dvdx_p + dvdx_m) / 2.0;
            vdvdy = v(ij) * (dvdy_p + dvdy_m) / 2.0;
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
void stepIntermediateSemiImplicitDense(Eigen::VectorXd& u_star, Eigen::VectorXd& v_star, const Eigen::VectorXd& u, const Eigen::VectorXd& v, const GridInfo& Grid, const BoundaryConditions& BC, const ProblemInformation& Problem)
{
    // indexing
    int ij;
    int ip1, im1; // (i+1,j) and (i-1,j)
    int jp1, jm1; // (i,j+1) and (i,j-1)
    // grid steps
    double dx_p, dx_m;
    double dy_p, dy_m;
    // laplacian coefficients
    double cx_p, cx_m, cx_0;
    double cy_p, cy_m, cy_0;
    // first order derivatives
    double dudx_p, dudx_m;
    double dudy_p, dudy_m;
    double dvdx_p, dvdx_m;
    double dvdy_p, dvdy_m;
    // nonlinear first order terms
    double ududx, vdudy;
    double udvdx, vdvdy;
    // system
    int NTOTAL = Grid.NX * Grid.NY; // TODO: not 3d yet lol
    Eigen::MatrixXd Au = Eigen::MatrixXd::Zero(NTOTAL, NTOTAL);
    Eigen::MatrixXd Av = Eigen::MatrixXd::Zero(NTOTAL, NTOTAL);
    Eigen::VectorXd bu = Eigen::VectorXd::Zero(NTOTAL);
    Eigen::VectorXd bv = Eigen::VectorXd::Zero(NTOTAL);

    for (int jdx = 0; jdx < Grid.NY; jdx++) {
        for (int idx = 0; idx < Grid.NX; idx++) {
            ij = jdx * Grid.NX + idx;
            im1 = jdx * Grid.NX + idx - 1;
            ip1 = jdx * Grid.NX + idx + 1;
            jm1 = (jdx - 1) * Grid.NX + idx;
            jp1 = (jdx + 1) * Grid.NX + idx;
            //
            if (BC.velocity_type[ij] == BoundaryConditionType::DIRECHLET) {
                Au(ij, ij) = 1.0;
                Av(ij, ij) = 1.0;
                bu(ij) = BC.u[ij];
                bv(ij) = BC.v[ij];
                continue;
            }
            // set grid steps
            dx_p = Grid.X(idx + 1) - Grid.X(idx);
            dx_m = Grid.X(idx) - Grid.X(idx - 1);
            dy_p = Grid.Y(jdx + 1) - Grid.Y(jdx);
            dy_m = Grid.Y(jdx) - Grid.Y(jdx - 1);
            // calculate laplacian coefficents
            cx_p = 2.0 / (dx_p * (dx_p + dx_m));
            cx_m = 2.0 / (dx_m * (dx_p + dx_m));
            cx_0 = -(cx_p + cx_m);
            cy_p = 2.0 / (dy_p * (dy_p + dy_m));
            cy_m = 2.0 / (dy_m * (dy_p + dy_m));
            cy_0 = -(cy_p + cy_m);
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
            ududx = u(ij) * (dudx_p + dudx_m) / 2.0;
            vdudy = v(ij) * (dudy_p + dudy_m) / 2.0;
            udvdx = u(ij) * (dvdx_p + dvdx_m) / 2.0;
            vdvdy = v(ij) * (dvdy_p + dvdy_m) / 2.0;
            // set values in coefficient matrix
            Au(ij, ij) = 1 / Problem.dt - (Problem.Properties.mu / Problem.Properties.rho) * (cx_0 + cy_0);
            Au(ij, ip1) = -(Problem.Properties.mu / Problem.Properties.rho) * cx_p;
            Au(ij, im1) = -(Problem.Properties.mu / Problem.Properties.rho) * cx_m;
            Au(ij, jp1) = -(Problem.Properties.mu / Problem.Properties.rho) * cy_p;
            Au(ij, jm1) = -(Problem.Properties.mu / Problem.Properties.rho) * cy_m;
            Av(ij, ij) = 1 / Problem.dt - (Problem.Properties.mu / Problem.Properties.rho) * (cx_0 + cy_0);
            Av(ij, ip1) = -(Problem.Properties.mu / Problem.Properties.rho) * cx_p;
            Av(ij, im1) = -(Problem.Properties.mu / Problem.Properties.rho) * cx_m;
            Av(ij, jp1) = -(Problem.Properties.mu / Problem.Properties.rho) * cy_p;
            Av(ij, jm1) = -(Problem.Properties.mu / Problem.Properties.rho) * cy_m;
            // set values in RHS
            bu(ij) = u(ij) / Problem.dt - ududx - vdudy;
            bv(ij) = v(ij) / Problem.dt - udvdx - vdvdy;
        }
    }

    // solve for intermediate velocity
    u_star = Au.partialPivLu().solve(bu);
    v_star = Av.partialPivLu().solve(bv);
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
    // grid steps
    double dx_p, dx_m;
    double dy_p, dy_m;
    // laplacian coefficients
    double cx_p, cx_m, cx_0;
    double cy_p, cy_m, cy_0;
    // first order derivatives
    double dudx_p, dudx_m;
    double dudy_p, dudy_m;
    double dvdx_p, dvdx_m;
    double dvdy_p, dvdy_m;
    // nonlinear first order terms
    double ududx, vdudy;
    double udvdx, vdvdy;
    // system
    int NTOTAL = Grid.NX * Grid.NY; // TODO: not 3d yet lol
    // Eigen::MatrixXd Au = Eigen::MatrixXd::Zero(NTOTAL, NTOTAL);
    // Eigen::MatrixXd Av = Eigen::MatrixXd::Zero(NTOTAL, NTOTAL);

    Eigen::SparseMatrix<double, Eigen::RowMajor> Au(NTOTAL, NTOTAL);
    std::vector<Eigen::Triplet<double>> utriplets;
    utriplets.reserve(5 * NTOTAL); // 2D problem has ~5 non-zero entries per row

    Eigen::SparseMatrix<double, Eigen::RowMajor> Av(NTOTAL, NTOTAL);
    std::vector<Eigen::Triplet<double>> vtriplets;
    vtriplets.reserve(5 * NTOTAL); // 2D problem has ~5 non-zero entries per row

    Eigen::VectorXd bu = Eigen::VectorXd::Zero(NTOTAL);
    Eigen::VectorXd bv = Eigen::VectorXd::Zero(NTOTAL);

    for (int jdx = 0; jdx < Grid.NY; jdx++) {
        for (int idx = 0; idx < Grid.NX; idx++) {
            ij = jdx * Grid.NX + idx;
            im1 = jdx * Grid.NX + idx - 1;
            ip1 = jdx * Grid.NX + idx + 1;
            jm1 = (jdx - 1) * Grid.NX + idx;
            jp1 = (jdx + 1) * Grid.NX + idx;
            //
            if (BC.velocity_type[ij] == BoundaryConditionType::DIRECHLET) {
                utriplets.push_back(Eigen::Triplet<double>(ij, ij, 1.0));
                vtriplets.push_back(Eigen::Triplet<double>(ij, ij, 1.0));
                bu(ij) = BC.u[ij];
                bv(ij) = BC.v[ij];
                continue;
            }
            // set grid steps
            dx_p = Grid.X(idx + 1) - Grid.X(idx);
            dx_m = Grid.X(idx) - Grid.X(idx - 1);
            dy_p = Grid.Y(jdx + 1) - Grid.Y(jdx);
            dy_m = Grid.Y(jdx) - Grid.Y(jdx - 1);
            // calculate laplacian coefficents
            cx_p = 2.0 / (dx_p * (dx_p + dx_m));
            cx_m = 2.0 / (dx_m * (dx_p + dx_m));
            cx_0 = -(cx_p + cx_m);
            cy_p = 2.0 / (dy_p * (dy_p + dy_m));
            cy_m = 2.0 / (dy_m * (dy_p + dy_m));
            cy_0 = -(cy_p + cy_m);
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
            ududx = u(ij) * (dudx_p + dudx_m) / 2.0;
            vdudy = v(ij) * (dudy_p + dudy_m) / 2.0;
            udvdx = u(ij) * (dvdx_p + dvdx_m) / 2.0;
            vdvdy = v(ij) * (dvdy_p + dvdy_m) / 2.0;
            // set values in coefficient matrix
            double au_ij = 1 / Problem.dt - (Problem.Properties.mu / Problem.Properties.rho) * (cx_0 + cy_0);
            double au_ip1 = -(Problem.Properties.mu / Problem.Properties.rho) * cx_p;
            double au_im1 = -(Problem.Properties.mu / Problem.Properties.rho) * cx_m;
            double au_jp1 = -(Problem.Properties.mu / Problem.Properties.rho) * cy_p;
            double au_jm1 = -(Problem.Properties.mu / Problem.Properties.rho) * cy_m;
            utriplets.push_back(Eigen::Triplet<double>(ij, ij, au_ij));
            utriplets.push_back(Eigen::Triplet<double>(ij, ip1, au_ip1));
            utriplets.push_back(Eigen::Triplet<double>(ij, im1, au_im1));
            utriplets.push_back(Eigen::Triplet<double>(ij, jp1, au_jp1));
            utriplets.push_back(Eigen::Triplet<double>(ij, jm1, au_jm1));
            //
            double av_ij = 1 / Problem.dt - (Problem.Properties.mu / Problem.Properties.rho) * (cx_0 + cy_0);
            double av_ip1 = -(Problem.Properties.mu / Problem.Properties.rho) * cx_p;
            double av_im1 = -(Problem.Properties.mu / Problem.Properties.rho) * cx_m;
            double av_jp1 = -(Problem.Properties.mu / Problem.Properties.rho) * cy_p;
            double av_jm1 = -(Problem.Properties.mu / Problem.Properties.rho) * cy_m;
            vtriplets.push_back(Eigen::Triplet<double>(ij, ij, av_ij));
            vtriplets.push_back(Eigen::Triplet<double>(ij, ip1, av_ip1));
            vtriplets.push_back(Eigen::Triplet<double>(ij, im1, av_im1));
            vtriplets.push_back(Eigen::Triplet<double>(ij, jp1, av_jp1));
            vtriplets.push_back(Eigen::Triplet<double>(ij, jm1, av_jm1));
            // set values in RHS
            bu(ij) = u(ij) / Problem.dt - ududx - vdudy;
            bv(ij) = v(ij) / Problem.dt - udvdx - vdvdy;
        }
    }
    // build coefficient matrix
    Au.setFromTriplets(utriplets.begin(), utriplets.end());
    Av.setFromTriplets(vtriplets.begin(), vtriplets.end());
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>> usolver;
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>> vsolver;
    // usolver.setMaxIterations(10000);
    // vsolver.setMaxIterations(10000);
    // usolver.setTolerance(1e-14);
    // vsolver.setTolerance(1e-14);
    usolver.compute(Au);
    vsolver.compute(Av);
    // solve for intermediate velocity
    u_star = usolver.solve(bu);
    v_star = vsolver.solve(bv);

    if (usolver.info() != Eigen::Success) {
        std::cerr << "Warning: Pressure solver did not converge. Iterations: "
                  << usolver.iterations() << ", Error: " << usolver.error() << std::endl;
    }
    if (vsolver.info() != Eigen::Success) {
        std::cerr << "Warning: Pressure solver did not converge. Iterations: "
                  << vsolver.iterations() << ", Error: " << vsolver.error() << std::endl;
    }
}