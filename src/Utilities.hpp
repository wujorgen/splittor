#ifndef UTILITIES
#define UTILITIES

#include <Eigen/Dense>

Eigen::VectorXd ConvertFieldToVector(const Eigen::ArrayXXd& field);

Eigen::ArrayXXd ConvertVectorToField(const Eigen::VectorXd& vector, const int& NX, const int& NY);

#endif