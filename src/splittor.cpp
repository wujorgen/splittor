#include <iostream>

#include "InputProcessor.hpp"
#include "InputReader.hpp"
#include "Types.hpp"

using namespace std;

int main()
{
    ProblemInformation Problem;
    GridInfo Grid;
    BoundaryConditions BC;
    FluidProperties Properties;

    int grid_reader_error = readGridFile(Grid, BC);
    if (grid_reader_error != 0) {
        std::cout << "A problem occured when reading the grid file. Exiting..." << std::endl;
        return grid_reader_error;
    }
    int problem_reader_error = readProblemInformation(Problem);
    if (problem_reader_error != 0) {
        std::cout << "A problem occured when reading the problem information file. Exiting..." << std::endl;
        return problem_reader_error;
    }

    return 0;
}