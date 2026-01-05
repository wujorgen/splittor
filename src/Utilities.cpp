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

/**
 * @brief Computes the max divergence of a flow field.
 *
 * @param u: vector store of flow field of size NX, NY
 * @param v: vector store of flow field of size NX, NY
 * @param Grid: grid information
 */
double computeMaxDivergence(const Eigen::VectorXd& u, const Eigen::VectorXd& v,
    const GridInfo& Grid, const bool& debug)
{
    double max_div = 0.0;
    double mean_div = 0.0;
    int count = 0;
    // indexing
    int ij;
    int ip1, im1; // (i+1,j) and (i-1,j)
    int jp1, jm1; // (i,j+1) and (i,j-1)
    // grid steps
    double dx_p, dx_m;
    double dy_p, dy_m;
    // first derivatives
    double dudx;
    double dvdy;
    double div;
    for (int jdx = 1; jdx < Grid.NY - 1; jdx++) {
        for (int idx = 1; idx < Grid.NX - 1; idx++) {
            ij = jdx * Grid.NX + idx;
            ip1 = jdx * Grid.NX + idx + 1;
            im1 = jdx * Grid.NX + idx - 1;
            jp1 = (jdx + 1) * Grid.NX + idx;
            jm1 = (jdx - 1) * Grid.NX + idx;

            dx_p = Grid.X(idx + 1) - Grid.X(idx);
            dx_m = Grid.X(idx) - Grid.X(idx - 1);
            dy_p = Grid.Y(jdx + 1) - Grid.Y(jdx);
            dy_m = Grid.Y(jdx) - Grid.Y(jdx - 1);

            dudx = 0.5 * ((u(ip1) - u(ij)) / dx_p + (u(ij) - u(im1)) / dx_m);
            dvdy = 0.5 * ((v(jp1) - v(ij)) / dy_p + (v(ij) - v(jm1)) / dy_m);

            div = dudx + dvdy;
            max_div = std::max(max_div, std::abs(div));
            mean_div += std::abs(div);
            count++;
            // std::cout << "idx: " << idx << ", jdx: " << jdx << ", div: " << div << std::endl;
        }
    }
    mean_div /= count;
    if (debug) {
        std::cout << "Max divergence: " << max_div
                  << ", Mean divergence: " << mean_div << std::endl;
    }
    return max_div;
}

/**
 * @brief Computes the max divergence of a flow field.
 *
 * @param u: vector store of flow field of size NX, NY
 * @param v: vector store of flow field of size NX, NY
 * @param Grid: grid information
 */
double computeMeanDivergence(const Eigen::VectorXd& u, const Eigen::VectorXd& v,
    const GridInfo& Grid, const bool& debug)
{
    double max_div = 0.0;
    double mean_div = 0.0;
    int count = 0;
    // indexing
    int ij;
    int ip1, im1; // (i+1,j) and (i-1,j)
    int jp1, jm1; // (i,j+1) and (i,j-1)
    // grid steps
    double dx_p, dx_m;
    double dy_p, dy_m;
    // first derivatives
    double dudx;
    double dvdy;
    double div;
    for (int jdx = 1; jdx < Grid.NY - 1; jdx++) {
        for (int idx = 1; idx < Grid.NX - 1; idx++) {
            ij = jdx * Grid.NX + idx;
            ip1 = jdx * Grid.NX + idx + 1;
            im1 = jdx * Grid.NX + idx - 1;
            jp1 = (jdx + 1) * Grid.NX + idx;
            jm1 = (jdx - 1) * Grid.NX + idx;

            dx_p = Grid.X(idx + 1) - Grid.X(idx);
            dx_m = Grid.X(idx) - Grid.X(idx - 1);
            dy_p = Grid.Y(jdx + 1) - Grid.Y(jdx);
            dy_m = Grid.Y(jdx) - Grid.Y(jdx - 1);

            dudx = 0.5 * ((u(ip1) - u(ij)) / dx_p + (u(ij) - u(im1)) / dx_m);
            dvdy = 0.5 * ((v(jp1) - v(ij)) / dy_p + (v(ij) - v(jm1)) / dy_m);

            div = dudx + dvdy;
            max_div = std::max(max_div, std::abs(div));
            mean_div += std::abs(div);
            count++;
            // std::cout << "idx: " << idx << ", jdx: " << jdx << ", div: " << div << std::endl;
        }
    }
    mean_div /= count;
    if (debug) {
        std::cout << "Max divergence: " << max_div
                  << ", Mean divergence: " << mean_div << std::endl;
    }
    return mean_div;
}

/**
 * @brief Computes divergence of a flow field.
 *
 * This version of the function attemps to take advantage of Eigen3's vectorization.
 */
DivergenceResult computeDivergence(const Eigen::VectorXd& u, const Eigen::VectorXd& v,
    const GridInfo& Grid, bool debug)
{
    // 1. Map flat vectors to 2D RowMajor views
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> U(u.data(), Grid.NY, Grid.NX);
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> V(v.data(), Grid.NY, Grid.NX);

    int innerNY = Grid.NY - 2;
    int innerNX = Grid.NX - 2;

    // 2. Mesh Spacing: delta = (x[i+1] - x[i-1])
    // We convert to .array() immediately to allow coefficient-wise math later
    Eigen::ArrayXd deltaX = (Grid.X.segment(2, innerNX) - Grid.X.segment(0, innerNX)).array();
    Eigen::ArrayXd deltaY = (Grid.Y.segment(2, innerNY) - Grid.Y.segment(0, innerNY)).array();

    // 3. Finite Difference numerator
    // We convert the Matrix blocks to Arrays to perform the division
    auto dudx_num = (U.block(1, 2, innerNY, innerNX) - U.block(1, 0, innerNY, innerNX)).array();
    auto dvdy_num = (V.block(2, 1, innerNY, innerNX) - V.block(0, 1, innerNY, innerNX)).array();

    // 4. Perform the non-uniform division (Broadcasting)
    // dudx: Divide each row by deltaX
    // dvdy: Divide each column by deltaY
    Eigen::ArrayXXd dudx = dudx_num.rowwise() / deltaX.transpose();
    Eigen::ArrayXXd dvdy = dvdy_num.colwise() / deltaY;

    // 5. Compute absolute divergence field
    Eigen::ArrayXXd div = (dudx + dvdy).abs();

    DivergenceResult res;
    res.max_div = div.maxCoeff();
    res.mean_div = div.mean();

    if (debug) {
        std::printf("Max Div: %.6e, Mean Div: %.6e\n", res.max_div, res.mean_div);
    }

    return res;
}