#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "Types.hpp"
#include "Utilities.hpp"

void writeToCSVFile(const Eigen::VectorXd& velocity, const GridInfo& Grid, const std::string fname)
{
    std::ofstream file(fname);
    file << std::setprecision(5);
    // indexing
    int ij;
    int ip1, im1; // (i+1,j) and (i-1,j)
    int jp1, jm1; // (i,j+1) and (i,j-1)
    //
    for (int jdx = 0; jdx < Grid.NY; jdx++) {
        for (int idx = 0; idx < Grid.NX; idx++) {
            // indexing
            ij = jdx * Grid.NX + idx;
            //
            file << velocity(ij);
            if (idx < Grid.NX - 1) {
                file << ", ";
            }
        }
        file << std::endl;
    }
}