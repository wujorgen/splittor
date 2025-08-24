#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xview.hpp>

using Eigen::ArrayXXd;
using Eigen::last;
using Eigen::seq;

void SolvePoissonGrid2D(xt::xarray<double>& phi, const xt::xarray<double>& b, const double dx, const double dy, const double ETOLP, const double MAXITERP) {
    /* Solves the Poisson equation (\nabla \cdot \nabla P = b). */
    // TODO: will need a way to specify boundary condition inputs
    double ErrorNorm;
    int NX = phi.shape(0);
    int NY = phi.shape(1);
    xt::xarray<double> phi_next = xt::zeros<double>({NX, NY});

    for (int itr_p = 0; itr_p < MAXITERP; itr_p++) {
        phi_next.fill(0.0);

        xt::view(phi_next, xt::range(1, -1), xt::range(1, -1)) =
            ((xt::view(phi, xt::range(2, NX), xt::range(1, -1)) + xt::view(phi, xt::range(0, -2), xt::range(1, -1))) * pow(dy, 2) +
             (xt::view(phi, xt::range(1, -1), xt::range(2, NY)) + xt::view(phi, xt::range(1, -1), xt::range(0, -2))) * pow(dx, 2) -
             xt::view(b, xt::range(1, -1), xt::range(1, -1)) * pow(dx, 2) * pow(dy, 2)) /
            (2 * (pow(dx, 2) + pow(dy, 2)));

        ErrorNorm = xt::linalg::norm(phi_next - phi);
        // std::cout << "PoissonGrid2D, Itr = " << itr_p << ", Error Norm: " << ErrorNorm << std::endl;
        std::swap(phi, phi_next);
        // phi = phi_next;

        if (ErrorNorm <= ETOLP) {
            break;
        }
    }
}

void SolvePoissonGrid2D_Eigen(ArrayXXd& phi, const ArrayXXd& b, const double dx, const double dy, const double ETOLP,
                              const double MAXITERP) {
    /* Solves the Poisson equation (\nabla \cdot \nabla P = b). */
    // TODO: will need a way to specify boundary condition inputs
    double ErrorNorm;
    ArrayXXd phi_next(phi.rows(), phi.cols());

    for (int itr_p = 0; itr_p < MAXITERP; itr_p++) {
        phi_next.setZero();

        phi_next(seq(1, last - 1), seq(1, last - 1)) =
            ((phi(seq(2, last), seq(1, last - 1)) + phi(seq(0, last - 2), seq(1, last - 1))) * pow(dy, 2) +
             (phi(seq(1, last - 1), seq(2, last)) + phi(seq(1, last - 1), seq(0, last - 2))) * pow(dx, 2) -
             b(seq(1, last - 1), seq(1, last - 1)) * pow(dx, 2) * pow(dy, 2)) /
            (2 * (pow(dx, 2) + pow(dy, 2)));

        // TODO: if there are pressure boundary conditions, enforce them here. They will come from a boundary condition object

        ErrorNorm = (phi_next - phi).matrix().norm();
        // std::cout << "PoissonGrid2D_Eigen, Itr = " << itr_p << ", Error Norm: " << ErrorNorm << std::endl;
        phi = phi_next;

        if (ErrorNorm <= ETOLP) {
            break;
        }
    }
}