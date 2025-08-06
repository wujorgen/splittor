#include <Eigen/Dense>
#include <cmath>

using Eigen::ArrayXXd;
using Eigen::last;
using Eigen::seq;

void SolvePoisson(ArrayXXd& phi, const ArrayXXd& b, const double dx, const double dy, const double ETOLP,
                  const double MAXITERP) {
    /* Solves the Poisson equation (\nabla \cdot \nabla P = b). */
    // TODO: ETOLP and MAXITERP should come from convergence property struct
    // TODO: dx and dy should come from a grid property struct
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
        phi = phi_next;

        if (ErrorNorm <= ETOLP) {
            break;
        }
    }
}