#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include "InputProcessor.hpp"
#include "Types.hpp"

/**
 * @brief Applies boundary conditions
 */
void applyBoundaryConditions2D(Eigen::VectorXd& u, Eigen::VectorXd& v, const BoundaryConditions& BC, const GridInfo& Grid)
{
    int ij;
    for (int jdx = 0; jdx < Grid.NY; jdx++) {
        for (int idx = 0; idx < Grid.NX; idx++) {
            ij = jdx * Grid.NX + idx;
            if (BC.type[ij] == BoundaryConditionType::VELOCITY_DIRECHLET) {
                u(ij) = BC.u[ij];
                v(ij) = BC.v[ij];
                continue;
            }
        }
    }
}

/**
 * @brief Calculates distances between grid points.
 *
 * TODO: for some reason this doesn't work right.
 * @param Grid
 */
void calcGridSteps(GridInfo& Grid)
{
    if (Grid.NX > 0)
        Grid.dx = Grid.X(Eigen::seq(1, Eigen::last)) - Grid.X(Eigen::seq(0, Eigen::last - 1));
    if (Grid.NY > 0)
        Grid.dy = Grid.Y(Eigen::seq(1, Eigen::last)) - Grid.Y(Eigen::seq(0, Eigen::last - 1));
    if (Grid.NZ > 0)
        Grid.dz = Grid.Z(Eigen::seq(1, Eigen::last)) - Grid.Z(Eigen::seq(0, Eigen::last - 1));
}

/**
 * @brief Checks if a grid is Cartesian.
 *
 * Loops over flattened lists of grid points, indexed by ijk, to check if a grid is Cartesian.
 * @param Grid
 * @param x
 * @param y
 * @param z
 * @param tolerance
 * @return true or false
 */
bool checkGridIsCartesian(const GridInfo& Grid, const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, double tolerance)
{
    bool isCartesian = true;
    // Check that x-coordinates only vary with i (constant along j and k)
    for (int k = 0; k < Grid.NZ; k++) {
        for (int j = 0; j < Grid.NY; j++) {
            for (int i = 0; i < Grid.NX; i++) {
                int ijk = k * (Grid.NX * Grid.NY) + j * Grid.NX + i;
                int ijk_ref = 0 * (Grid.NX * Grid.NY) + 0 * Grid.NX + i; // same i, j=0, k=0

                if (std::abs(x[ijk] - x[ijk_ref]) > tolerance) {
                    std::cerr << "Grid is not Cartesian: x varies with j or k" << std::endl;
                    isCartesian = false;
                    break;
                }
            }
        }
    }
    // Check that y-coordinates only vary with j (constant along i and k)
    for (int k = 0; k < Grid.NZ; k++) {
        for (int j = 0; j < Grid.NY; j++) {
            for (int i = 0; i < Grid.NX; i++) {
                int ijk = k * (Grid.NX * Grid.NY) + j * Grid.NX + i;
                int ijk_ref = 0 * (Grid.NX * Grid.NY) + j * Grid.NX + 0; // same j, i=0, k=0

                if (std::abs(y[ijk] - y[ijk_ref]) > tolerance) {
                    std::cerr << "Grid is not Cartesian: y varies with i or k" << std::endl;
                    isCartesian = false;
                    break;
                }
            }
        }
    }
    // Check that z-coordinates only vary with k (constant along i and j)
    for (int k = 0; k < Grid.NZ; k++) {
        for (int j = 0; j < Grid.NY; j++) {
            for (int i = 0; i < Grid.NX; i++) {
                int ijk = k * (Grid.NX * Grid.NY) + j * Grid.NX + i;
                int ijk_ref = k * (Grid.NX * Grid.NY) + 0 * Grid.NX + 0; // same k, i=0, j=0

                if (std::abs(z[ijk] - z[ijk_ref]) > tolerance) {
                    std::cerr << "Grid is not Cartesian: z varies with i or j" << std::endl;
                    std::cerr << "\t" << z[ijk] << " " << z[ijk_ref] << std::endl;
                    isCartesian = false;
                    break;
                }
            }
        }
    }
    return isCartesian;
}