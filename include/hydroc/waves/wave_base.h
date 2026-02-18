/*********************************************************************
 * @file  wave_base.h
 * @brief Base classes shared by all wave models.
 *********************************************************************/

#ifndef HYDROC_WAVES_WAVE_BASE_H
#define HYDROC_WAVES_WAVE_BASE_H

#include <Eigen/Dense>

//// RADU - why is this not in some namespace?
////      - should these classes be DLL-exported (e.g., to allow AddWave(NoWave))?
////        if not, why is this a public headers?

enum class WaveMode {
    noWaveCIC = 0,
    regular   = 1,
    irregular = 2
};

/**
 * @brief Abstract base class for all wave models.
 *
 * Coordinate conventions:
 *   - Waves propagate in the +X direction (unidirectional).
 *   - Z is vertical (positive upward), with z = mwl_ at mean water level.
 *   - Y is horizontal, perpendicular to wave propagation.
 *
 * Units:
 *   - Positions: meters [m]
 *   - Time: seconds [s]
 *   - Elevation η: meters [m]
 *   - Gradients ∂η/∂x, ∂η/∂y: dimensionless [m/m]
 *   - Velocities: meters per second [m/s]
 *   - Accelerations: meters per second squared [m/s²]
 */
class WaveBase {
  public:
    virtual ~WaveBase() = default;

    virtual void Initialize()                                                                               = 0;
    virtual Eigen::VectorXd GetForceAtTime(double t)                                                        = 0;
    virtual WaveMode GetWaveMode()                                                                          = 0;
    virtual double GetElevation(const Eigen::Vector3d& position, double time)                               = 0;
    virtual Eigen::Vector3d GetVelocity(const Eigen::Vector3d& position, double time, double elevation)     = 0;
    virtual Eigen::Vector3d GetAcceleration(const Eigen::Vector3d& position, double time, double elevation) = 0;

    Eigen::Vector3d GetVelocity(const Eigen::Vector3d& position, double time);
    Eigen::Vector3d GetAcceleration(const Eigen::Vector3d& position, double time);

    void SetNumBodies(unsigned int n) { num_bodies_ = n; }
    unsigned int GetNumBodies() const { return num_bodies_; }

    double mwl_           = 0.0;
    double g_             = 9.81;
    double water_depth_   = 0.0;
    bool wave_stretching_ = true;

  protected:
    unsigned int num_bodies_ = 0;
};

class NoWave : public WaveBase {
  public:
    NoWave() = default;

    void Initialize() override {}
    Eigen::VectorXd GetForceAtTime(double t) override;
    WaveMode GetWaveMode() override { return mode_; }
    double GetElevation(const Eigen::Vector3d&, double) override { return 0.0; }
    Eigen::Vector3d GetVelocity(const Eigen::Vector3d&, double, double) override { return Eigen::Vector3d(0.0, 0.0, 0.0); }
    Eigen::Vector3d GetAcceleration(const Eigen::Vector3d&, double, double) override { return Eigen::Vector3d(0.0, 0.0, 0.0); }

  private:
    const WaveMode mode_ = WaveMode::noWaveCIC;
};

#endif  // HYDROC_WAVES_WAVE_BASE_H

