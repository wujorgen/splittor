#include <Eigen/Dense>

#include "CalcPressure.hpp"
#include "Types.hpp"

void calcPressure(Eigen::VectorXd& p_np1, const Eigen::VectorXd& u_star, const Eigen::VectorXd& v_star, const GridInfo& Grid, const ProblemInformation& Problem)
{
    // indexing
    int ij;
    int ip1, im1; // (i+1,j) and (i-1,j)
    int jp1, jm1; // (i,j+1) and (i,j-1)

    // coefficient matrix and RHS
    Eigen::MatrixXd A(Grid.NX * Grid.NY, Grid.NX * Grid.NY);
    Eigen::VectorXd b(Grid.NX * Grid.NY);
    A.setZero();
    b.setZero();

    for (int jdx = 0; jdx < Grid.NY; jdx++) {
        for (int idx = 0; idx < Grid.NX; idx++) {
            ij = jdx * Grid.NX + idx;
            im1 = jdx * Grid.NX + idx - 1;
            ip1 = jdx * Grid.NX + idx + 1;
            jm1 = (jdx - 1) * Grid.NX + idx;
            jp1 = (jdx + 1) * Grid.NX + idx;
            // handle boundaries
            // grid steps
            // first order terms
            // build coefficient matrix
            // build RHS
        }
    }
    // solve
    // TODO: gotta check. https://libeigen.gitlab.io/eigen/docs-nightly/group__TutorialLinearAlgebra.html
    p_np1 = A.partialPivLu().solve(b);
}