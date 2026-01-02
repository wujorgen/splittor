#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <iomanip> // needed for setprecision
#include <iostream>
#include <vector>

#include "../src/CalcIntermediateVelocity.hpp"
#include "../src/CalcPressure.hpp"
#include "../src/InputReader.hpp"
#include "../src/PISO.hpp"
#include "../src/Types.hpp"
#include "../src/Utilities.hpp"
#include "../src/VelocityProjection.hpp"

TEST(OperatorSplitting, Chorin)
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
    Eigen::VectorXd u_next = Eigen::VectorXd::Zero(Grid.NX * Grid.NY);
    Eigen::VectorXd v_next = Eigen::VectorXd::Zero(Grid.NX * Grid.NY);
    Eigen::VectorXd p_star = Eigen::VectorXd::Zero(Grid.NX * Grid.NY);

    applyVelocityBoundaryConditions2D(u, v, BC, Grid);

    // below is a Chorin projection step

    // solve momentum for intermediate velocity
    // stepIntermediateExplicit(u_star, v_star, u, v, Grid, BC, Problem);
    stepIntermediateSemiImplicit(u_star, v_star, u, v, Grid, BC, Problem);

    // solve poisson equation for pressure field
    calcPressure(p_star, u_star, v_star, Grid, BC, Problem);

    // projection step
    projectVelocity(u_next, v_next, u_star, v_star, p_star, Grid, BC, Problem);

    // printing and checking
    bool print = false;
    Eigen::MatrixXd u_next_mat = ConvertVectorToField(u_next, Grid.NX, Grid.NY);
    for (int jdx = 0; jdx < Grid.NY; jdx++) {
        for (int idx = 0; idx < Grid.NX; idx++) {
            if (print) {
                std::cout << std::fixed << std::setprecision(5) << u_next_mat(idx, jdx) << ", ";
            }
            switch (Grid.NY - jdx) {
            case 1:
                EXPECT_NEAR(u_next_mat(idx, jdx), 2, 1e-6);
                break;
            // case 2:
            //     if (idx != 0 && idx != Grid.NX - 1) {
            //         EXPECT_NEAR(u_star_mat(idx, jdx), 0.2, 1e-6);
            //     }
            //     break;
            default:
                break;
            }
        }
        if (print) {
            std::cout << std::endl;
        }
    }
}

TEST(OperatorSplitting, PISO_STEP)
{
    GridInfo Grid;
    BoundaryConditions BC;
    ProblemInformation Problem;
    int grid_reader_error = readGridFile(Grid, BC, "testdata/lid_grid.txt");
    int problem_reader_error = readProblemInformation(Problem, "testdata/lid_probleminformation.txt");

    Eigen::VectorXd u = Eigen::VectorXd::Zero(Grid.NX * Grid.NY);
    Eigen::VectorXd v = Eigen::VectorXd::Zero(Grid.NX * Grid.NY);
    Eigen::VectorXd p = Eigen::VectorXd::Zero(Grid.NX * Grid.NY);

    Eigen::VectorXd u_s1 = Eigen::VectorXd::Zero(Grid.NX * Grid.NY);
    Eigen::VectorXd v_s1 = Eigen::VectorXd::Zero(Grid.NX * Grid.NY);
    Eigen::VectorXd p_s1 = Eigen::VectorXd::Zero(Grid.NX * Grid.NY);
    Eigen::VectorXd u_s2 = Eigen::VectorXd::Zero(Grid.NX * Grid.NY);
    Eigen::VectorXd v_s2 = Eigen::VectorXd::Zero(Grid.NX * Grid.NY);

    Eigen::VectorXd u_next = Eigen::VectorXd::Zero(Grid.NX * Grid.NY);
    Eigen::VectorXd v_next = Eigen::VectorXd::Zero(Grid.NX * Grid.NY);
    Eigen::VectorXd p_next = Eigen::VectorXd::Zero(Grid.NX * Grid.NY);

    applyVelocityBoundaryConditions2D(u, v, BC, Grid);

    stepPISO(u_next, v_next, p_next, u, v, p, Grid, BC, Problem);

    // printing and checking
    bool print = true;
    Eigen::MatrixXd u_next_mat = ConvertVectorToField(u_next, Grid.NX, Grid.NY);
    for (int jdx = 0; jdx < Grid.NY; jdx++) {
        for (int idx = 0; idx < Grid.NX; idx++) {
            if (print) {
                std::cout << std::fixed << std::setprecision(5) << u_next_mat(idx, jdx) << ", ";
            }
            switch (Grid.NY - jdx) {
            case 1:
                // EXPECT_NEAR(u_next_mat(idx, jdx), 2, 1e-6);
                break;
            // case 2:
            //     if (idx != 0 && idx != Grid.NX - 1) {
            //         EXPECT_NEAR(u_star_mat(idx, jdx), 0.2, 1e-6);
            //     }
            //     break;
            default:
                break;
            }
        }
        if (print) {
            std::cout << std::endl;
        }
    }
}