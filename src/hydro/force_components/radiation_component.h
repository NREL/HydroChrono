/*********************************************************************
 * @file  radiation_component.h
 * @brief Radiation damping force component (velocity history convolution).
 *********************************************************************/

#ifndef HYDRO_FORCE_COMPONENTS_RADIATION_COMPONENT_H
#define HYDRO_FORCE_COMPONENTS_RADIATION_COMPONENT_H

#include <hydroc/core/force_component.h>
#include <hydroc/core/system_state.h>
#include <hydroc/radiation/radiation_types.h>
#include "../radiation/radiation_rirf_convolution.h"
#include "../radiation/radiation_rirf_processing.h"
#include <Eigen/Dense>
#include <vector>
#include <memory>

class HydroData;

namespace hydrochrono::hydro {

/**
 * @brief Radiation damping force component (RIRF convolution).
 * 
 * Computes radiation damping forces via convolution of RIRF kernels
 * with body velocity history. Supports optional kernel processing
 * (smoothing and/or tapering) via RadiationKernelProcessing.
 */
class RadiationComponent : public IHydroForceComponent {
public:
    /**
     * @brief Constructor.
     * 
     * @param file_info Reference to HydroData containing RIRF kernels
     * @param num_bodies Number of bodies in the system
     * @param rirf_steps Number of time steps in RIRF (truncated to [0, T] by
     *                    HydroSystem::CreateRadiationComponent() when truncation is set)
     * @param rirf_time_vector Time vector for RIRF (same truncation applies)
     * @param rirf_width_vector Width vector for RIRF (for trapezoidal integration)
     * @param kernel_processing Composable smoothing/tapering options (defaults are no-ops)
     * @param diagnostics_output_dir Directory for diagnostics output (CSV files)
     */
    RadiationComponent(const HydroData& file_info,
                  int num_bodies,
                  int rirf_steps,
                  const Eigen::VectorXd& rirf_time_vector,
                  const Eigen::VectorXd& rirf_width_vector,
                  const RadiationKernelProcessing& kernel_processing,
                  const std::string& diagnostics_output_dir);

    HydroComponentType Type() const override { return HydroComponentType::Radiation; }

    void Compute(const SystemState& state,
                double time,
                BodyForces& inout_forces) override;

private:
    const HydroData& file_info_;
    int num_bodies_;
    RadiationKernelProcessing kernel_processing_;
    std::string diagnostics_output_dir_;
    
    std::unique_ptr<RadiationRirfConvolution> convolution_;
    Eigen::VectorXd rirf_time_vector_;
    Eigen::VectorXd rirf_width_vector_;
    
    bool rirf_processed_ready_ = false;
    std::vector<Eigen::Tensor<double, 3>> rirf_processed_;
    
    static constexpr int kDofPerBody = 6;
    
    void EnsureProcessedRIRF();
    double GetRIRFval(int row, int col, int st);
};

}  // namespace hydrochrono::hydro

#endif  // HYDRO_FORCE_COMPONENTS_RADIATION_COMPONENT_H

