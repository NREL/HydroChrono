/*********************************************************************
 * @file  wave_base.cpp
 * @brief Base wave model implementations.
 *********************************************************************/

#include <hydroc/waves/wave_base.h>

Eigen::Vector3d WaveBase::GetVelocity(const Eigen::Vector3d& position, double time) {
    // Elevation needed only if stretching=true
    double elevation = wave_stretching_ ? GetElevation(position, time) : 0.0;
    return GetVelocity(position, time, elevation);
}

Eigen::Vector3d WaveBase::GetAcceleration(const Eigen::Vector3d& position, double time) {
    // Elevation needed only if stretching=true
    double elevation = wave_stretching_ ? GetElevation(position, time) : 0.0;
    return GetAcceleration(position, time, elevation);
}

NoWave::NoWave() : num_bodies_(1) {}

NoWave::NoWave(unsigned int num_b) : num_bodies_(num_b) {}

NoWave::NoWave(const NoWaveParams& params) : num_bodies_(params.num_bodies_) {}

Eigen::VectorXd NoWave::GetForceAtTime(double) {
    unsigned int dof = num_bodies_ * 6;
    Eigen::VectorXd f(dof);
    for (int i = 0; i < f.size(); i++) {
        f[i] = 0.0;
    }
    return f;
}
