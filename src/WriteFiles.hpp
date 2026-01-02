#ifndef WRITE_FILES
#define WRITE_FILES

#include <Eigen/Dense>
#include "Types.hpp"

void writeToCSVFile(const Eigen::VectorXd& velocity, const GridInfo& Grid, const std::string fname);

#endif