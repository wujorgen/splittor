#include <iostream>
#include <Eigen/Dense>

#include "Types.hpp"

void CalcGridSteps(GridInfo& Grid) {
    if (Grid.NX > 0)
        Grid.dx = Grid.X(Eigen::seq(1, Eigen::last)) - Grid.X(Eigen::seq(0, Eigen::last - 1));
    if (Grid.NY > 0)
        Grid.dy = Grid.Y(Eigen::seq(1, Eigen::last)) - Grid.Y(Eigen::seq(0, Eigen::last - 1));
    if (Grid.NZ > 0)
        Grid.dz = Grid.Z(Eigen::seq(1, Eigen::last)) - Grid.Z(Eigen::seq(0, Eigen::last - 1));
}