/*********************************************************************
 * @file  wave_base.h
 * @brief Base classes shared by all wave models.
 *********************************************************************/

#ifndef HYDROC_WAVES_WAVE_BASE_H
#define HYDROC_WAVES_WAVE_BASE_H

#include <Eigen/Dense>

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

    virtual void Initialize()                                   = 0;
    virtual Eigen::VectorXd GetForceAtTime(double t)            = 0;
    virtual WaveMode         GetWaveMode()                      = 0;
    virtual double           GetElevation(const Eigen::Vector3d& position, double time)    = 0;
    virtual Eigen::Vector3d  GetVelocity(const Eigen::Vector3d& position, double time)     = 0;
    virtual Eigen::Vector3d  GetAcceleration(const Eigen::Vector3d& position, double time) = 0;

    /**
     * @brief Returns the free-surface gradient (∂η/∂x, ∂η/∂y) at the given position and time.
     *
     * Used for computing surface normals in visualization. The normal vector can be
     * constructed as: n = normalize(−∂η/∂x, −∂η/∂y, 1).
     *
     * @param position World coordinates (x, y, z) in meters. Only x, y are used.
     * @param time     Simulation time in seconds.
     * @return Eigen::Vector2d containing (∂η/∂x, ∂η/∂y), dimensionless.
     *
     * @note Default implementation returns (0, 0) for flat surface (NoWave).
     * @note Current wave models are unidirectional (+X), so ∂η/∂y = 0.
     */
    virtual Eigen::Vector2d GetElevationGradientXY(const Eigen::Vector3d& position, double time) const {
        (void)position;
        (void)time;
        return Eigen::Vector2d::Zero();
    }

    double mwl_         = 0.0;
    double g_           = 9.81;
    double water_depth_ = 0.0;
};

class NoWave : public WaveBase {
  public:
    NoWave();
    explicit NoWave(unsigned int num_b);

    void Initialize() override {}
    Eigen::VectorXd GetForceAtTime(double t) override;
    WaveMode GetWaveMode() override { return mode_; }
    double GetElevation(const Eigen::Vector3d&, double) override { return 0.0; }
    Eigen::Vector3d GetVelocity(const Eigen::Vector3d&, double) override { return Eigen::Vector3d(0.0, 0.0, 0.0); }
    Eigen::Vector3d GetAcceleration(const Eigen::Vector3d&, double) override { return Eigen::Vector3d(0.0, 0.0, 0.0); }

  private:
    unsigned int num_bodies_ = 0;
    const WaveMode mode_ = WaveMode::noWaveCIC;
};

#endif  // HYDROC_WAVES_WAVE_BASE_H

