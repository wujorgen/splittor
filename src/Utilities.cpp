#include <Eigen/Dense>
// #include <unsupported/Eigen/CXX11/Tensor>

#include "Utilities.hpp"

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
 * @brief
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
