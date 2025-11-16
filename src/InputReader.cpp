#include "InputReader.hpp"

#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "InputProcessor.hpp"
#include "InputReader.hpp"
#include "Types.hpp"

int ReadInputFile(ProblemInformation& Problem, GridInfo& Grid, BoundaryConditions& BC, const std::string& fname) {
    // TODO: There must be a better way to do this. Until then - this is what we have.
    std::ifstream file(fname);  // TODO: that's not the right filename.

    if (!file.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return -1;
    }
    bool HEADER_FOUND = false;
    std::string line;

    std::string str_a;
    std::string str_b;
    std::string str_c;
    std::string str_d;
    std::string str_e;
    int i_a, i_b, i_c, i_d, i_e;
    double d_a, d_b, d_c, d_d, d_e;

    // Header must be first line of grid input.
    // Formatted as "NX 3 NY 2 NZ 1 D 2"
    std::getline(file, line);
    std::stringstream ss(line);
    if (ss >> str_a >> i_a >> str_b >> i_b >> str_c >> i_c >> str_d >> i_d) {
        Grid.NX = i_a;
        Grid.NY = i_b;
        Grid.NZ = i_c;
        // std::cout << i_a << "," << i_b << "," << i_c << "," << i_d << std::endl;
        HEADER_FOUND = true;
        std::cout << "SEPERATION" << std::endl;
    } else {
        std::cerr << "Header not detected." << std::endl;
        return -1;
    }

    int i = 0;
    int j = 0;
    int k = 0;
    int ijk = 0;
    std::vector<double> x_locs;
    std::vector<double> y_locs;
    std::vector<double> z_locs;
    std::vector<double> bc_types;
    std::vector<double> bc_vals;
    while (std::getline(file, line)) {
        // Skip empty lines
        if (line.empty()) {
            continue;
        }
        std::stringstream ss(line);
        if (ss >> d_a >> d_b >> d_c >> i_d >> d_e) {
            // Try to read 3 doubles, 1 int, and another double
            // std::cout << "Grid point: " << d_a << " " << d_b << " " << d_c << std::endl;
            // std::cout << "\tBC Type: " << i_d << " BC Value: " << d_e << std::endl;
            // keep track of processed grid points using ijk indexing
            ijk = k * (Grid.NX * Grid.NY) + j * Grid.NX + i;
            //
            // std::cout << "\ti:" << i << " j:" << j << " k:" << k << " ijk: " << ijk << std::endl;
            x_locs.push_back(d_a);
            y_locs.push_back(d_b);
            z_locs.push_back(d_c);
            bc_types.push_back(i_d);
            bc_vals.push_back(d_e);
            //
            ijk++;  // Increment the flat index and decompose for next valid line
            i = ijk % Grid.NX;
            j = (ijk / Grid.NX) % Grid.NY;
            k = ijk / (Grid.NX * Grid.NY);
        } else {
            // Otherwise try to read as a line of strings
            // Reset the stringstream to read from beginning
            // These lines wil probably just be ignored as comments
            ss.clear();
            ss.seekg(0);
            std::vector<std::string> words;
            std::string word;
            while (ss >> word) {
                words.push_back(word);
            }
            std::cout << "Line of strings: ";
            for (const auto& w : words) {
                std::cout << w << " ";
            }
            std::cout << std::endl;
        }
    }
    file.close();
    for (int index = 0; index < ijk; index++) {
        std::cout << x_locs[index] << " " << y_locs[index] << " " << z_locs[index] << std::endl;
    }

    Eigen::VectorXd eigen_vec = Eigen::Map<Eigen::VectorXd>(bc_types.data(), bc_types.size());  // TODO: syntax for mapping vector to eigen vector (not eigenvector)
    if (!GridIsCartesian(Grid, x_locs, y_locs, z_locs)) {
        std::cerr << "The input grid is not cartesian." << std::endl;
        return 1;
    }
    return 0;
}

//int main() {
//    ProblemInformation Problem;
//    GridInfo Grid;
//    BoundaryConditions BC;
//    FluidProperties Properties;
//    ReadInputFile(Problem, Grid, BC);
//    std::cout << Grid.NX << " x " << Grid.NY << " x " << Grid.NZ << std::endl;
//}