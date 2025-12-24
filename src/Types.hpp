#ifndef TYPES
#define TYPES

#include <string>
#include <Eigen/Dense>

// Store fluid properties
struct FluidProperties {
    double mu = 1;
    double rho = 1;
};

// Store convergence params
struct ConvergenceSettings {
    double relax = 1;
};

// Store problem information
struct ProblemInformation {
    int Dimensions = 2;
    double dt = 0.01;
    double EndTime = 10;
    std::string Mode = "steady";
    FluidProperties Properties;
    ConvergenceSettings Convergence;
};

// Store a flow field
struct FlowField {
    Eigen::VectorXd u;
    Eigen::VectorXd v;
    Eigen::VectorXd w;
};

// Store grid information
struct GridInfo {
    //
    int NX;
    int NY;
    int NZ;
    //
    Eigen::VectorXd X;
    Eigen::VectorXd Y;
    Eigen::VectorXd Z;
    //
    Eigen::VectorXd dx;
    Eigen::VectorXd dy;
    Eigen::VectorXd dz;
};

// Store boundary conditions
struct BoundaryConditions {
    // store if point is boundary point
    // 0 = no, 1 = velocity direchlet, 2 = pressure direchlet, 3 = velocity neumann, 4 = pressure direchlet
    Eigen::VectorXi XBC;
    Eigen::VectorXi YBC;
    Eigen::VectorXi ZBC;
    Eigen::VectorXi PBC;
    // boundary condition value
    Eigen::VectorXd u;
    Eigen::VectorXd v;
    Eigen::VectorXd w;
    Eigen::VectorXd p;
};

#endif