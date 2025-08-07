#include <gtest/gtest.h>

#include <Eigen/Dense>

#include "../src/SolverPoisson.hpp"

TEST(SolverPoisson, Grid2D) {
    int NX = 11;
    int NY = 11;
    double XMIN = 0;
    double XMAX = 2;
    double YMIN = 0;
    double YMAX = 2;

    double dx = (XMAX - XMIN) / (NX - 1);
    double dy = (YMAX - YMIN) / (NY - 1);

    Eigen::ArrayXXd phi(NX, NY);
    Eigen::ArrayXXd b(NX, NY);
    phi.setZero();
    b.setZero();

    // TODO: initialize correct solution
    Eigen::ArrayXXd phiCorrect(NX, NY);

    // set source terms
    b(int(NX/4), int(NX/4)) = 100;
    b(int(3*NY/4), int(3*NY/4)) = -100;

    SolvePoissonGrid2D(phi, b, dx, dy);

}