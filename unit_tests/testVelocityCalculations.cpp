#include <gtest/gtest.h>

#include <Eigen/Dense>

#include "../src/CalcIntermediateVelocity.hpp"
#include "../src/InputReader.hpp"
#include "../src/Types.hpp"
#include "../src/Utilities.hpp"

TEST(VelocityCalculation, ExplicitStep)
{
    GridInfo Grid;
    BoundaryConditions BC;
    ProblemInformation Problem;
    int grid_reader_error = readGridFile(Grid, BC, "testdata/lid_grid.txt");
    int problem_reader_error = readProblemInformation(Problem, "testdata/lid_probleminformation.txt");
    Eigen::VectorXd u = Eigen::VectorXd::Zero(Grid.NX * Grid.NY);
    Eigen::VectorXd v = Eigen::VectorXd::Zero(Grid.NX * Grid.NY);
    Eigen::VectorXd u_star = Eigen::VectorXd::Zero(Grid.NX * Grid.NY);
    Eigen::VectorXd v_star = Eigen::VectorXd::Zero(Grid.NX * Grid.NY);

    // TODO: need a seperate routine to apply boundary conditions to u and v
    int ij;
    for (int jdx = 0; jdx < Grid.NY; jdx++) {
        for (int idx = 0; idx < Grid.NX; idx++) {
            ij = jdx * Grid.NX + idx;
            if (jdx == Grid.NY - 1) {
                u(ij) = 2.0;
                v(ij) = 0.0;
                continue;
            } else if (jdx == 0 or idx == 0 or idx == Grid.NX - 1) {
                u(ij) = 0.0;
                v(ij) = 0.0;
                continue;
            }
        }
    }

    stepIntermediateExplicit(u_star, v_star, u, v, Grid, BC, Problem);

    Eigen::MatrixXd u_star_mat = ConvertVectorToField(u_star, Grid.NX, Grid.NY);
    for (int jdx = 0; jdx < Grid.NY; jdx++) {
        for (int idx = 0; idx < Grid.NX; idx++) {
            // std::cout << u_star_mat(idx, jdx) << ", ";
            switch (Grid.NY - jdx) {
            case 1:
                EXPECT_NEAR(u_star_mat(idx, jdx), 2, 1e-6);
                break;
            case 2:
                if (idx != 0 && idx != Grid.NX - 1) {
                    EXPECT_NEAR(u_star_mat(idx, jdx), 0.2, 1e-6);
                }
                break;
            default:
                break;
            }
        }
        // std::cout << std::endl;
    }
    for (int ij = 0; ij < Grid.NX * Grid.NY; ij++) {
        EXPECT_NEAR(v_star(ij), 0.0, 1e-6);
    }
}

TEST(VelocityCalculation, SemiImplicitStep)
{
    EXPECT_TRUE(false);
}