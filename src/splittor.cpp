#include <fstream>
#include <iostream>

#include "Constants.hpp"
#include "InputProcessor.hpp"
#include "InputReader.hpp"
#include "SolveSteadyState.hpp"
#include "SolveTransient.hpp"
#include "Types.hpp"
#include "Utilities.hpp"
#include "WriteFiles.hpp"

using namespace std;

int main()
{
    std::cout << "The number of threads that Eigen3 can see is: " << Eigen::nbThreads() << std::endl;

    ProblemInformation Problem;
    GridInfo Grid;
    BoundaryConditions BC;
    FluidProperties Properties;

    int grid_reader_error = readGridFile(Grid, BC, Constants::FNAME_GRID);
    if (grid_reader_error != 0) {
        std::cerr << "A problem occured when reading the grid file. Exiting..." << std::endl;
        return grid_reader_error;
    }

    int problem_reader_error = readProblemInformation(Problem, Constants::FNAME_PROBLEM_INFORMATION);
    if (problem_reader_error != 0) {
        std::cerr << "A problem occured when reading the problem information file. Exiting..." << std::endl;
        return problem_reader_error;
    }

    // open and clear all output files? - TODO this should be its own routine
    std::ofstream(FNAME_OUTPUT_U);
    std::ofstream(FNAME_OUTPUT_V);
    // std::ofstream(FNAME_OUTPUT_W);
    std::ofstream(FNAME_OUTPUT_P);
    std::ofstream(FNAME_LOG);

    int NTOTAL = Grid.NX * Grid.NY * Grid.NZ;
    Eigen::VectorXd u = Eigen::VectorXd::Zero(NTOTAL);
    Eigen::VectorXd v = Eigen::VectorXd::Zero(NTOTAL);
    Eigen::VectorXd w = Eigen::VectorXd::Zero(NTOTAL);
    Eigen::VectorXd p = Eigen::VectorXd::Zero(NTOTAL);

    if (Problem.Dimensions == 2) {
        applyVelocityBoundaryConditions2D(u, v, BC, Grid);
    }

    if (Problem.Mode == std::string("steady")) {
        solveSteadyStateProblem(u, v, p, Grid, BC, Problem);
        writeToCSVFile(u, Grid, Constants::FNAME_OUTPUT_U);
        writeToCSVFile(v, Grid, Constants::FNAME_OUTPUT_V);
        writeToCSVFile(p, Grid, Constants::FNAME_OUTPUT_P);
    } else if (Problem.Mode == std::string("transient")) {
    } else {
        std::cerr << "An invalid problem type has been requested. Exiting..." << std::endl;
        return 1;
    }

    return 0;
}