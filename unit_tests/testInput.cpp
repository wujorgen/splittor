#include <gtest/gtest.h>

#include <Eigen/Dense>

#include "../src/InputReader.hpp"
#include "../src/InputProcessor.hpp"
#include "../src/Types.hpp"

TEST(Input, Reader) {
    ProblemInformation Problem;
    GridInfo Grid;
    BoundaryConditions BC;
    FluidProperties Properties;
    ReadInputFile(Problem, Grid, BC, "testdata/grid.txt");
    std::cout << Grid.NX << " x " << Grid.NY << " x " << Grid.NZ << std::endl;
}

TEST(Input, Processor) {}

TEST(InputProcessor, CalcGridSteps) {
    GridInfo Grid;
    Grid.NX = 11;
    Grid.NY = 21;
    Grid.NZ = 5;

    Grid.X = Eigen::VectorXd::LinSpaced(Grid.NX, 0, 1);
    Grid.Y = Eigen::VectorXd::LinSpaced(Grid.NY, 0, 1);

    Eigen::VectorXd zvec = Eigen::VectorXd::Zero(Grid.NZ);
    zvec << 0, 0.1, 0.3, 0.6, 1.0;
    Grid.Z = zvec;

    Eigen::VectorXd dzvec = Eigen::VectorXd::Zero(Grid.NZ - 1);
    dzvec << 0.1, 0.2, 0.3, 0.4;

    CalcGridSteps(Grid);

    for (int i = 0; i < Grid.NX - 1; i++) {
        EXPECT_NEAR(Grid.dx(i), 0.1, 1e-6);
    }
    for (int j = 0; j < Grid.NY - 1; j++) {
        EXPECT_NEAR(Grid.dy(j), 0.05, 1e-6);
    }
    for (int k = 0; k < Grid.NZ - 1; k++) {
        EXPECT_NEAR(Grid.dz(k), dzvec(k), 1e-6);
    }
}