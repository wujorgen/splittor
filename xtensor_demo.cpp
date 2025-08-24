#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xio.hpp>

// g++ xtensor_demo.cpp -lblas -llapack
// Ubuntu/apt does not need explicit include call for xtensor
// on Mac, when xtensor and xtensor-blas are installed via Brew, explicit linking might be needed.

int main() {
    xt::xarray<double> A = {{1.0, 0.0}, {0.0, 1.0}};
    xt::xarray<double> b = {1.0, 1.0};  // 1D is fine
    xt::xarray<double> x = xt::linalg::solve(A, b);
    std::cout << "A: " << A << std::endl;
    std::cout << "b: " << b << std::endl;
    std::cout << "Solution to Ax=b: " << x << std::endl;

    xt::xarray<double> matrix = {{1, 2}, {3, 4}};  // 2x2
    xt::xarray<double> v = {1, 2};            // 1D

    auto result1 = xt::linalg::dot(matrix, v);     // v as column: A×v = 2×1
    auto result2 = xt::linalg::dot(v, matrix);     // v as row: v×A = 1×2

    // xtensor does not auto promote - both of the above are just shape (2,)
    std::cout << "matrix times v: " << result1 << " with shape " << result1.shape()[0] << " by " << result1.shape()[1] << std::endl;
    std::cout << "v times matrix: " << result2 << " with shape " << result2.shape()[0] << " by " << result2.shape()[1] << std::endl;

    return 0;
}