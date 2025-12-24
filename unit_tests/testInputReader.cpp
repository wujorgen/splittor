#include <gtest/gtest.h>

#include <Eigen/Dense>

#include "../src/InputProcessor.hpp"
#include "../src/InputReader.hpp"
#include "../src/Types.hpp"

TEST(InputReader, Grid332)
{
    ProblemInformation Problem;
    GridInfo Grid;
    BoundaryConditions BC;
    // FluidProperties Properties;
    int InputReaderError = ReadGridFile(Problem, Grid, BC, "testdata/grid332.txt");
    if (false)
    { // to view grid points during unit testing
        std::cout << "Grid Points X: ";
        for (int i = 0; i < Grid.NX; i++)
        {
            std::cout << Grid.X[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "Grid Points Y: ";
        for (int j = 0; j < Grid.NY; j++)
        {
            std::cout << Grid.Y[j] << ", ";
        }
        std::cout << std::endl;
        std::cout << "Grid Points Z: ";
        for (int k = 0; k < Grid.NZ; k++)
        {
            std::cout << Grid.Z[k] << ", ";
        }
        std::cout << std::endl;
    }
    if (false)
    { // to view grid deltas during unit testing
        CalcGridSteps(Grid);
        std::cout << "Grid Deltas X: ";
        for (int i = 0; i < Grid.NX - 1; i++)
        {
            std::cout << Grid.dx[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "Grid Deltas Y: ";
        for (int j = 0; j < Grid.NY - 1; j++)
        {
            std::cout << Grid.dy[j] << ", ";
        }
        std::cout << std::endl;
        std::cout << "Grid Deltas Z: ";
        for (int k = 0; k < Grid.NZ - 1; k++)
        {
            std::cout << Grid.dz[k] << ", ";
        }
        std::cout << std::endl;
    }
    EXPECT_EQ(InputReaderError, 0);
    // std::cout << Grid.NX << " x " << Grid.NY << " x " << Grid.NZ << std::endl;
}

TEST(InputReader, GridNoHeader)
{
    ProblemInformation Problem;
    GridInfo Grid;
    BoundaryConditions BC;
    // FluidProperties Properties;
    int InputReaderError = ReadGridFile(Problem, Grid, BC, "testdata/gridnoheader.txt");
    EXPECT_EQ(InputReaderError, 1);
}

TEST(InputReader, ProblemInfo) {}