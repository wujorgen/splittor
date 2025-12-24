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

int ReadGridFile(ProblemInformation &Problem, GridInfo &Grid, BoundaryConditions &BC, const std::string &fname)
{
    // TODO: There must be a better way to do this. Until then - this is what we have.
    std::ifstream file(fname); // TODO: that's not the right filename.

    if (!file.is_open())
    {
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
    int i_nx, i_ny, i_nz, i_nd, i_bctype;
    double d_xval, d_yval, d_zval, d_uval, d_vval, d_pval;

    // Header must be first line of grid input.
    // Formatted as "NX 3 NY 2 NZ 1 D 2"
    std::getline(file, line);
    std::stringstream ss(line);
    if (ss >> str_a >> i_nx >> str_b >> i_ny >> str_c >> i_nz >> str_d >> i_nd)
    {
        Grid.NX = i_nx;
        Grid.NY = i_ny;
        Grid.NZ = i_nz;
        // std::cout << i_a << "," << i_b << "," << i_c << "," << i_d << std::endl;
        HEADER_FOUND = true;
        // std::cout << "SEPERATION" << std::endl;
    }
    else
    {
        std::cerr << "Header not detected." << std::endl;
        return 1;
    }

    int i = 0;
    int j = 0;
    int k = 0;
    int ijk = 0;
    std::vector<double> x_locs;
    std::vector<double> y_locs;
    std::vector<double> z_locs;
    std::vector<double> bc_types;
    std::vector<double> u_vals;
    std::vector<double> v_vals;
    std::vector<double> p_vals;
    while (std::getline(file, line))
    {
        // Skip empty lines
        if (line.empty())
        {
            continue;
        }
        std::stringstream ss(line);
        // Formatting of rest of file:
        // X  Y  Z  BCTYPE  UVAL  VVAL  PVAL
        // d  d  d  i       d     d     d
        if (ss >> d_xval >> d_yval >> d_zval >> i_bctype >> d_uval >> d_vval >> d_pval)
        {
            // std::cout << "Grid point: " << d_xval << " " << d_yval << " " << d_zval << std::endl;
            // std::cout << "\tBC Type: " << i_d << " BC Value: " << d_uval << "," << d_vval << "," << d_pval << std::endl;
            //  keep track of processed grid points using ijk indexing
            ijk = k * (Grid.NX * Grid.NY) + j * Grid.NX + i;
            //
            // std::cout << "\ti:" << i << " j:" << j << " k:" << k << " ijk: " << ijk << std::endl;
            x_locs.push_back(d_xval);
            y_locs.push_back(d_yval);
            z_locs.push_back(d_zval);
            bc_types.push_back(i_bctype);
            u_vals.push_back(d_uval);
            v_vals.push_back(d_vval);
            p_vals.push_back(d_pval);
            //
            ijk++; // Increment the flat index and decompose for next valid line
            i = ijk % Grid.NX;
            j = (ijk / Grid.NX) % Grid.NY;
            k = ijk / (Grid.NX * Grid.NY);
        }
        else
        {
            // Otherwise try to read as a line of strings
            // Reset the stringstream to read from beginning
            // These lines wil probably just be ignored as comments
            ss.clear();
            ss.seekg(0);
            std::vector<std::string> words;
            std::string word;
            while (ss >> word)
            {
                words.push_back(word);
            }
            // std::cout << "Line of strings: ";
            // for (const auto& w : words) {
            //     std::cout << w << " ";
            // }
            // std::cout << std::endl;
        }
    }

    file.close();

    // for (int index = 0; index < ijk; index++) {
    //     std::cout << x_locs[index] << " " << y_locs[index] << " " << z_locs[index] << std::endl;
    // }

    if (!GridIsCartesian(Grid, x_locs, y_locs, z_locs))
    {
        std::cerr << "The input grid is not cartesian." << std::endl;
        return 2;
    }

    // if the grid has been confirmed as cartesian, get the unique grid points for x, y, and z
    std::vector<double> unique_x;
    std::vector<double> unique_y;
    std::vector<double> unique_z;
    for (int i = 0; i < Grid.NX; i++)
    {
        ijk = 0 * (Grid.NX * Grid.NY) + 0 * Grid.NX + i;
        unique_x.push_back(x_locs[ijk]);
    }
    for (int j = 0; j < Grid.NX; j++)
    {
        ijk = 0 * (Grid.NX * Grid.NY) + j * Grid.NX + 0;
        unique_y.push_back(y_locs[ijk]);
    }
    for (int k = 0; k < Grid.NX; k++)
    {
        ijk = k * (Grid.NX * Grid.NY) + 0 * Grid.NX + 0;
        unique_z.push_back(z_locs[ijk]);
    }
    Grid.X = Eigen::Map<Eigen::VectorXd>(unique_x.data(), unique_x.size()); // map std::vector to VectorXd
    Grid.Y = Eigen::Map<Eigen::VectorXd>(unique_y.data(), unique_y.size());
    Grid.Z = Eigen::Map<Eigen::VectorXd>(unique_z.data(), unique_z.size());

    return 0;
}

int ReadProblemInformation(ProblemInformation &Problem, GridInfo &Grid, BoundaryConditions &BC, const std::string &fname)
{
    return 0;
}

// int main() {
//     ProblemInformation Problem;
//     GridInfo Grid;
//     BoundaryConditions BC;
//     FluidProperties Properties;
//     ReadInputFile(Problem, Grid, BC);
//     std::cout << Grid.NX << " x " << Grid.NY << " x " << Grid.NZ << std::endl;
// }