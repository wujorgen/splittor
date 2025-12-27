#include <Eigen/Dense>
// #include <unsupported/Eigen/CXX11/Tensor>
#include <iostream>

#include "Types.hpp"
#include "Utilities.hpp"

/**
 * @brief Applies boundary conditions
 */
void applyVelocityBoundaryConditions2D(Eigen::VectorXd& u, Eigen::VectorXd& v, const BoundaryConditions& BC, const GridInfo& Grid)
{
    int ij;
    for (int jdx = 0; jdx < Grid.NY; jdx++) {
        for (int idx = 0; idx < Grid.NX; idx++) {
            ij = jdx * Grid.NX + idx;
            if (BC.velocity_type[ij] == BoundaryConditionType::DIRECHLET) {
                // std::cout << "setting velocity direchlet in aVBC2D " << ij << std::endl;
                u(ij) = BC.u[ij];
                v(ij) = BC.v[ij];
            }
        }
    }
}

/**
 * @brief Flattens a field to a vector. Iteration is done row first.
 *
 * @param field shape (i,j)
 * @return vector of shape (i*j)
 */
Eigen::VectorXd ConvertFieldToVector(const Eigen::ArrayXXd& field)
{
    int NX = field.rows();
    int NY = field.cols();
    Eigen::VectorXd vector(NX * NY);
    int ij = 0;
    for (int jdx = 0; jdx < NY; jdx++) {
        for (int idx = 0; idx < NX; idx++) {
            ij = jdx * NX + idx;
            vector(ij) = field(idx, jdx);
        }
    }
    return vector;
}

/**
 * @brief Unflattens a vector to a field. Iteration is done row first.
 *
 * @param vector shape (i*j)
 * @return matrix of shape (i,j)
 */
Eigen::ArrayXXd ConvertVectorToField(const Eigen::VectorXd& vector, const int& NX, const int& NY)
{
    Eigen::ArrayXXd field(NX, NY);
    int ij = 0;
    for (int jdx = 0; jdx < NY; jdx++) {
        for (int idx = 0; idx < NX; idx++) {
            ij = jdx * NX + idx;
            field(idx, jdx) = vector(ij);
        }
    }
    return field;
}
