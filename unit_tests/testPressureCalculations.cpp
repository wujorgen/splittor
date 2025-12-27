#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <iomanip>

#include "../src/CalcIntermediateVelocity.hpp"
#include "../src/CalcPressure.hpp"
#include "../src/InputReader.hpp"
#include "../src/Types.hpp"
#include "../src/Utilities.hpp"

TEST(PressureCalculation, PressureCalculation)
{
    GridInfo Grid;
    BoundaryConditions BC;
    ProblemInformation Problem;
    int grid_reader_error = readGridFile(Grid, BC, "testdata/lid_grid.txt");
    int problem_reader_error = readProblemInformation(Problem, "testdata/lid_probleminformation.txt");

    int NTOTAL = Grid.NX * Grid.NY;
    Eigen::VectorXd u = Eigen::VectorXd::Zero(NTOTAL);
    Eigen::VectorXd v = Eigen::VectorXd::Zero(NTOTAL);
    Eigen::VectorXd u_star = Eigen::VectorXd::Zero(NTOTAL);
    Eigen::VectorXd v_star = Eigen::VectorXd::Zero(NTOTAL);
    Eigen::VectorXd p_star = Eigen::VectorXd::Zero(NTOTAL);

    Eigen::MatrixXd p_star_check(Grid.NX, Grid.NY);
    p_star_check << -0.036, -0.036, -0.040, -0.049, -0.065, -0.090, -0.128, -0.190, -0.295, -0.441, 0.000,
        -0.036, -0.036, -0.040, -0.049, -0.065, -0.090, -0.128, -0.190, -0.295, -0.441, 0.000,
        -0.031, -0.031, -0.035, -0.043, -0.056, -0.076, -0.104, -0.142, -0.185, -0.192, 0.000,
        -0.023, -0.023, -0.026, -0.032, -0.041, -0.054, -0.070, -0.089, -0.100, -0.082, 0.000,
        -0.012, -0.012, -0.014, -0.017, -0.021, -0.027, -0.035, -0.042, -0.043, -0.031, 0.000,
        0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
        0.012, 0.012, 0.014, 0.017, 0.021, 0.027, 0.035, 0.042, 0.043, 0.031, 0.000,
        0.023, 0.023, 0.026, 0.032, 0.041, 0.054, 0.070, 0.089, 0.100, 0.082, 0.000,
        0.031, 0.031, 0.035, 0.043, 0.056, 0.076, 0.104, 0.142, 0.185, 0.192, 0.000,
        0.036, 0.036, 0.040, 0.049, 0.065, 0.090, 0.128, 0.190, 0.295, 0.441, 0.000,
        0.036, 0.036, 0.040, 0.049, 0.065, 0.090, 0.128, 0.190, 0.295, 0.441, 0.000;

    applyVelocityBoundaryConditions2D(u, v, BC, Grid);

    stepIntermediateSemiImplicit(u_star, v_star, u, v, Grid, BC, Problem);

    calcPressure(p_star, u_star, v_star, Grid, BC, Problem);

    Eigen::MatrixXd p_mat = ConvertVectorToField(p_star, Grid.NX, Grid.NY);
    for (int jdx = 0; jdx < Grid.NY; jdx++) {
        for (int idx = 0; idx < Grid.NX; idx++) {
            // std::cout << std::fixed << std::setprecision(3) << p_star_check(idx, jdx) << ", ";
            EXPECT_NEAR(p_mat(idx, jdx), p_star_check(idx, jdx), 1e-3);
        }
        // std::cout << std::endl;
    }
}