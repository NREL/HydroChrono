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

class WaveBase {
  public:
    virtual ~WaveBase() = default;

    virtual void Initialize()                                   = 0;
    virtual Eigen::VectorXd GetForceAtTime(double t)            = 0;
    virtual WaveMode         GetWaveMode()                      = 0;
    virtual double           GetElevation(const Eigen::Vector3d& position, double time)    = 0;
    virtual Eigen::Vector3d  GetVelocity(const Eigen::Vector3d& position, double time)     = 0;
    virtual Eigen::Vector3d  GetAcceleration(const Eigen::Vector3d& position, double time) = 0;

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

