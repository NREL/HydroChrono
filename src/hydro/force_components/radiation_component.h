/*********************************************************************
 * @file  radiation_component.h
 * @brief Radiation damping force component (velocity history convolution).
 *********************************************************************/

#ifndef HYDRO_FORCE_COMPONENTS_RADIATION_COMPONENT_H
#define HYDRO_FORCE_COMPONENTS_RADIATION_COMPONENT_H

#include <hydroc/core/force_component.h>
#include <hydroc/core/system_state.h>
#include "../radiation/radiation_rirf_convolution.h"
#include "../radiation/radiation_rirf_processing.h"
#include <Eigen/Dense>
#include <vector>
#include <memory>

class HydroData;

namespace hydrochrono::hydro {

/**
 * @brief Convolution mode for radiation damping.
 */
enum class RadiationConvolutionMode {
    Baseline,
    TaperedDirect
};

/**
 * @brief Radiation damping force component (RIRF convolution).
 * 
 * Computes radiation damping forces via convolution of RIRF kernels
 * with body velocity history. Supports Baseline and TaperedDirect modes.
 */
class RadiationComponent : public IHydroForceComponent {
public:
    /**
     * @brief Constructor.
     * 
     * @param file_info Reference to HydroData containing RIRF kernels
     * @param num_bodies Number of bodies in the system
     * @param rirf_steps Number of time steps in RIRF
     * @param rirf_time_vector Time vector for RIRF
     * @param rirf_width_vector Width vector for RIRF (for trapezoidal integration)
     * @param convolution_mode Convolution mode (Baseline or TaperedDirect)
     * @param tapered_opts Options for TaperedDirect preprocessing
     * @param diagnostics_output_dir Directory for diagnostics output (CSV files)
     */
    RadiationComponent(const HydroData& file_info,
                  int num_bodies,
                  int rirf_steps,
                  const Eigen::VectorXd& rirf_time_vector,
                  const Eigen::VectorXd& rirf_width_vector,
                  RadiationConvolutionMode convolution_mode,
                  const TaperedDirectOptions& tapered_opts,
                  const std::string& diagnostics_output_dir);

    /**
     * @brief Get the component type.
     * @return HydroComponentType::Radiation
     */
    HydroComponentType Type() const override { return HydroComponentType::Radiation; }

    /**
     * @brief Compute radiation damping force contribution.
     * 
     * Records velocities from state, performs RIRF convolution,
     * and adds radiation damping forces to inout_forces per body.
     * 
     * @param state Current system state (velocities needed for history)
     * @param time Current simulation time
     * @param inout_forces Force vector to add contribution to (one GeneralizedForce per body)
     */
    void Compute(const SystemState& state,
                double time,
                BodyForces& inout_forces) override;

private:
    const HydroData& file_info_;
    int num_bodies_;
    RadiationConvolutionMode convolution_mode_;
    TaperedDirectOptions tapered_opts_;
    std::string diagnostics_output_dir_;
    
    std::unique_ptr<RadiationRirfConvolution> convolution_;
    Eigen::VectorXd rirf_time_vector_;
    Eigen::VectorXd rirf_width_vector_;
    
    // Processed RIRF tensors (for TaperedDirect mode)
    bool rirf_processed_ready_;
    std::vector<Eigen::Tensor<double, 3>> rirf_processed_;
    
    static constexpr int kDofPerBody = 6;
    
    void EnsureProcessedRIRF();
    double GetRIRFval(int row, int col, int st);
};

}  // namespace hydrochrono::hydro

#endif  // HYDRO_FORCE_COMPONENTS_RADIATION_COMPONENT_H

