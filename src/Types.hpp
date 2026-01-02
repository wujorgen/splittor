#ifndef TYPES
#define TYPES

#include <Eigen/Dense>
#include <string>
#include <vector>

// Store fluid properties
struct FluidProperties {
    double mu = 1;
    double rho = 1;
};

// Store PISO convergence settings
struct PISOSettings {
    int MAX_NUM_ITERATIONS = 500;
    int PRINT_EVERY_N_ITERATIONS = 50;
    double TOL_DIVERGENCE = 1e-6;
    double TOL_U_REL = 1e-4;
    double TOL_V_REL = 1e-4;
    double TOL_W_REL = 1e-4;
};

struct SteadyStateSettings {
    double TOL_U_REL = 1e-4;
    double TOL_V_REL = 1e-4;
    double TOL_W_REL = 1e-4;
    double MAX_STEPS = 1000; // 1000
};

struct TransientSettings {
    double TOL_U_REL = 1e-4;
    double TOL_V_REL = 1e-4;
    double TOL_W_REL = 1e-4;
};

// Store all convergence params
struct ConvergenceSettings {
    double relax = 1;
    PISOSettings PISO;
    SteadyStateSettings Steady;
    TransientSettings Transient;
};

// Store problem information
struct ProblemInformation {
    int Dimensions = 2;
    double dt = 0.01;
    double EndTime = 10;
    int nt;
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

// 0 = no, 1 = direchlet, 2 = neumann
enum BoundaryConditionType {
    NONE,
    DIRECHLET,
    PRESSURE
};

// Store boundary conditions
struct BoundaryConditions {
    // store if point is boundary point
    std::vector<BoundaryConditionType> velocity_type;
    std::vector<BoundaryConditionType> pressure_type;
    // boundary condition values
    std::vector<double> u;
    std::vector<double> v;
    std::vector<double> w;
    // TODO: since i'm lazy - edges of domain all default to no pressure gradient.
    std::vector<double> p;
};

#endif