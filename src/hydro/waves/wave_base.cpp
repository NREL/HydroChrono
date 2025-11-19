/*********************************************************************
 * @file  wave_base.cpp
 * @brief Base wave model implementations.
 *********************************************************************/

#include "wave_base.h"

NoWave::NoWave() : num_bodies_(1) {}

NoWave::NoWave(unsigned int num_b) : num_bodies_(num_b) {}

Eigen::VectorXd NoWave::GetForceAtTime(double) {
    unsigned int dof = num_bodies_ * 6;
    Eigen::VectorXd f(dof);
    for (int i = 0; i < f.size(); i++) {
        f[i] = 0.0;
    }
    return f;
}

