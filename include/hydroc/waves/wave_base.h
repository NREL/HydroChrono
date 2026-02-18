/*********************************************************************
 * @file  wave_base.h
 * @brief Base classes shared by all wave models.
 *********************************************************************/

#ifndef HYDROC_WAVES_WAVE_BASE_H
#define HYDROC_WAVES_WAVE_BASE_H

#include <Eigen/Core>

// TODO: Wrap wave classes in a namespace (e.g., hydrochrono::waves).
// TODO: Determine if these classes need DLL-export for external wave model support.

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
    static constexpr int kDofsPerBody = 6;

    virtual ~WaveBase() = default;

    WaveBase(const WaveBase&)            = delete;
    WaveBase& operator=(const WaveBase&) = delete;

    virtual void Initialize()                                                                                     = 0;
    virtual Eigen::VectorXd GetForceAtTime(double t) const                                                        = 0;
    virtual WaveMode GetWaveMode() const                                                                          = 0;
    virtual double GetElevation(const Eigen::Vector3d& position, double time) const                               = 0;
    virtual Eigen::Vector3d GetVelocity(const Eigen::Vector3d& position, double time, double elevation) const     = 0;
    virtual Eigen::Vector3d GetAcceleration(const Eigen::Vector3d& position, double time, double elevation) const = 0;

    Eigen::Vector3d GetVelocity(const Eigen::Vector3d& position, double time) const;
    Eigen::Vector3d GetAcceleration(const Eigen::Vector3d& position, double time) const;

    void SetNumBodies(unsigned int n) { num_bodies_ = n; }
    unsigned int GetNumBodies() const { return num_bodies_; }

    double GetMWL() const { return mwl_; }
    double GetGravity() const { return g_; }
    double GetWaterDepth() const { return water_depth_; }
    bool GetWaveStretching() const { return wave_stretching_; }

  protected:
    WaveBase() = default;

    double mwl_           = 0.0;
    double g_             = 9.81;
    double water_depth_   = 0.0;
    bool wave_stretching_ = true;
    unsigned int num_bodies_ = 0;
};

class NoWave : public WaveBase {
  public:
    NoWave() = default;

    void Initialize() override {}
    Eigen::VectorXd GetForceAtTime(double t) const override;
    WaveMode GetWaveMode() const override { return mode_; }
    double GetElevation(const Eigen::Vector3d&, double) const override { return 0.0; }
    Eigen::Vector3d GetVelocity(const Eigen::Vector3d&, double, double) const override { return Eigen::Vector3d(0.0, 0.0, 0.0); }
    Eigen::Vector3d GetAcceleration(const Eigen::Vector3d&, double, double) const override { return Eigen::Vector3d(0.0, 0.0, 0.0); }

  private:
    static constexpr WaveMode mode_ = WaveMode::noWaveCIC;
};

#endif  // HYDROC_WAVES_WAVE_BASE_H

