#include <gtest/gtest.h>

#include <Eigen/Dense>

#include "../src/InputProcessor.hpp"
#include "../src/InputReader.hpp"
#include "../src/Types.hpp"

TEST(Input, Grid332) {
    ProblemInformation Problem;
    GridInfo Grid;
    BoundaryConditions BC;
    // FluidProperties Properties;
    int InputReaderError = ReadGridFile(Problem, Grid, BC, "testdata/grid332.txt");
    if (false) { // to view grid points during unit testing
        std::cout << "Grid Points X: ";
        for (int i = 0; i < Grid.NX; i++) {
            std::cout << Grid.X[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "Grid Points Y: ";
        for (int j = 0; j < Grid.NY; j++) {
            std::cout << Grid.Y[j] << ", ";
        }
        std::cout << std::endl;
        std::cout << "Grid Points Z: ";
        for (int k = 0; k < Grid.NZ; k++) {
            std::cout << Grid.Z[k] << ", ";
        }
        std::cout << std::endl;
    }
    if (false) { // to view grid deltas during unit testing
        CalcGridSteps(Grid);
        std::cout << "Grid Deltas X: ";
        for (int i = 0; i < Grid.NX - 1; i++) {
            std::cout << Grid.dx[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "Grid Deltas Y: ";
        for (int j = 0; j < Grid.NY - 1; j++) {
            std::cout << Grid.dy[j] << ", ";
        }
        std::cout << std::endl;
        std::cout << "Grid Deltas Z: ";
        for (int k = 0; k < Grid.NZ - 1; k++) {
            std::cout << Grid.dz[k] << ", ";
        }
        std::cout << std::endl;
    }
    EXPECT_EQ(InputReaderError, 0);
    // std::cout << Grid.NX << " x " << Grid.NY << " x " << Grid.NZ << std::endl;
}

TEST(Input, GridNoHeader) {
    ProblemInformation Problem;
    GridInfo Grid;
    BoundaryConditions BC;
    // FluidProperties Properties;
    int InputReaderError = ReadGridFile(Problem, Grid, BC, "testdata/gridnoheader.txt");
    EXPECT_EQ(InputReaderError, 1);
}

TEST(Input, GridNotCartesian) {
    // TODO
    EXPECT_TRUE(true);
}

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