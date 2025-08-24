#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xview.hpp>

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

    xt::xarray<double> phi = xt::zeros<double>({NX, NY});
    xt::xarray<double> b = xt::zeros<double>({NX, NY});

    xt::xarray<double> phiCorrect = {
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, -0.290308, -0.580616, -0.379162, -0.210571, -0.10912, -0.051637, -0.0202656, -0.00506552, 0, 0},
        {0, -0.580616, -1.65304, -0.725462, -0.354053, -0.174271, -0.0771886, -0.02436, 0, 0.00506552, 0},
        {0, -0.379162, -0.725462, -0.515684, -0.305906, -0.156787, -0.0584867, 0, 0.02436, 0.0202656, 0},
        {0, -0.210571, -0.354053, -0.305906, -0.197196, -0.0884858, 0, 0.0584867, 0.0771886, 0.051637, 0},
        {0, -0.10912, -0.174271, -0.156787, -0.0884858, 0, 0.0884858, 0.156787, 0.174271, 0.10912, 0},
        {0, -0.051637, -0.0771886, -0.0584867, 0, 0.0884858, 0.197196, 0.305906, 0.354053, 0.210571, 0},
        {0, -0.0202656, -0.02436, 0, 0.0584867, 0.156787, 0.305906, 0.515684, 0.725462, 0.379162, 0},
        {0, -0.00506552, 0, 0.02436, 0.0771886, 0.174271, 0.354053, 0.725462, 1.65304, 0.580616, 0},
        {0, 0, 0.00506552, 0.0202656, 0.051637, 0.10912, 0.210571, 0.379162, 0.580616, 0.290308, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};

    // set source terms
    b(int(NX / 4), int(NX / 4)) = 100;
    b(int(3 * NY / 4), int(3 * NY / 4)) = -100;

    SolvePoissonGrid2D(phi, b, dx, dy, 1e-4, 100);
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            // EXPECT_DOUBLE_EQ(phi(i,j), phiCorrect(i,j));
            EXPECT_NEAR(phi(i, j), phiCorrect(i, j), 1e-6);
        }
    }
}

TEST(SolverPoisson, Grid2D_Eigen) {
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

    phiCorrect << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, -0.290308, -0.580616, -0.379162, -0.210571, -0.10912, -0.051637, -0.0202656, -0.00506552, 0, 0,
        0, -0.580616, -1.65304, -0.725462, -0.354053, -0.174271, -0.0771886, -0.02436, 0, 0.00506552, 0,
        0, -0.379162, -0.725462, -0.515684, -0.305906, -0.156787, -0.0584867, 0, 0.02436, 0.0202656, 0,
        0, -0.210571, -0.354053, -0.305906, -0.197196, -0.0884858, 0, 0.0584867, 0.0771886, 0.051637, 0,
        0, -0.10912, -0.174271, -0.156787, -0.0884858, 0, 0.0884858, 0.156787, 0.174271, 0.10912, 0,
        0, -0.051637, -0.0771886, -0.0584867, 0, 0.0884858, 0.197196, 0.305906, 0.354053, 0.210571, 0,
        0, -0.0202656, -0.02436, 0, 0.0584867, 0.156787, 0.305906, 0.515684, 0.725462, 0.379162, 0,
        0, -0.00506552, 0, 0.02436, 0.0771886, 0.174271, 0.354053, 0.725462, 1.65304, 0.580616, 0,
        0, 0, 0.00506552, 0.0202656, 0.051637, 0.10912, 0.210571, 0.379162, 0.580616, 0.290308, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

    // set source terms
    b(int(NX / 4), int(NX / 4)) = 100;
    b(int(3 * NY / 4), int(3 * NY / 4)) = -100;

    SolvePoissonGrid2D_Eigen(phi, b, dx, dy, 1e-4, 100);
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            // EXPECT_DOUBLE_EQ(phi(i,j), phiCorrect(i,j));
            EXPECT_NEAR(phi(i, j), phiCorrect(i, j), 1e-6);
        }
    }

    /*
    The following python code confirms the above result. It was run using https://www.online-python.com/

    import numpy as np

    NX = 11
    NY = 11

    XMIN = 0
    XMAX = 2
    YMIN = 0
    YMAX = 2

    dx = (XMAX-XMIN)/(NX-1)
    dy = (YMAX-YMIN)/(NY-1)

    phi = np.zeros((NX,NY))
    b = np.zeros((NX,NY))
    b[NX//4, NX//4] = 100
    b[3*NY//4, 3*NY//4] = -100

    phi_next = np.zeros((NX, NY))
    for itr in range(0, 100):
        phi_next = np.zeros((NX, NY))
        phi_next[1:-1,1:-1] = ((phi[2:,1:-1] + phi[0:-2,1:-1]) * dy**2 +
                                (phi[1:-1,2:] + phi[1:-1,0:-2]) * dx**2 -
                                b[1:-1,1:-1]*(dx**2)*(dy**2)) / (2*(dx**2 + dy**2))
        errornorm = np.linalg.norm(phi_next-phi)
        print(f"Iteration: {itr}, {errornorm=}")
        phi = phi_next
        if errornorm <= 1e-4:
            break
    print(np.round(phi,3))
    */
}