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

    double mwl_           = 0.0;
    double g_             = 9.81;
    double water_depth_   = 0.0;
    bool wave_stretching_ = true;
};

struct NoWaveParams {
    unsigned int num_bodies_;
};

class NoWave : public WaveBase {
  public:
    NoWave();
    explicit NoWave(unsigned int num_b);
    explicit NoWave(const NoWaveParams& params);

    void Initialize() override {}
    Eigen::VectorXd GetForceAtTime(double t) override;
    WaveMode GetWaveMode() override { return mode_; }
    double GetElevation(const Eigen::Vector3d&, double) override { return 0.0; }
    Eigen::Vector3d GetVelocity(const Eigen::Vector3d&, double, double) override { return Eigen::Vector3d(0.0, 0.0, 0.0); }
    Eigen::Vector3d GetAcceleration(const Eigen::Vector3d&, double, double) override { return Eigen::Vector3d(0.0, 0.0, 0.0); }

  private:
    unsigned int num_bodies_;
    const WaveMode mode_ = WaveMode::noWaveCIC;
};

#endif  // HYDROC_WAVES_WAVE_BASE_H

