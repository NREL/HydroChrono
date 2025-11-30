/*********************************************************************
 * @file  radiation_ss_component.h
 * @brief Radiation damping force component (state-space approximation).
 *
 * This component provides an alternative to convolution-based radiation
 * damping. Instead of convolving the RIRF kernel with velocity history
 * (O(N) per timestep), it uses a state-space model with exponential and
 * oscillatory modes (O(1) per timestep).
 *
 * INITIALIZATION:
 *   - Reads RIRF kernels K_ij(t) from HydroData
 *   - Fits each kernel using RadiationStateSpaceFitter
 *   - Creates RadiationStateSpaceModel for each DOF pair
 *
 * RUNTIME:
 *   - Steps all models forward with current velocities
 *   - Accumulates radiation forces from all models
 *   - O(1) per timestep vs O(N) for convolution
 *
 * DOF COUPLING:
 *   For a multi-DOF system with N bodies (6N total DOFs), the full RIRF
 *   kernel K is a (6N×6N×T) tensor. Each K_ij(t) represents the force
 *   in DOF i due to unit velocity in DOF j.
 *
 *   This implementation creates a separate state-space model for each
 *   non-zero (i,j) pair. Most pairs are non-zero for single-body systems,
 *   but cross-body coupling depends on the BEM analysis.
 *
 * TYPICAL SPEEDUP:
 *   - 10-20x faster than convolution for simulations > 1000 steps
 *   - Speedup increases with simulation length
 *********************************************************************/

#ifndef HYDRO_FORCE_COMPONENTS_RADIATION_SS_COMPONENT_H
#define HYDRO_FORCE_COMPONENTS_RADIATION_SS_COMPONENT_H

#include <hydroc/core/force_component.h>
#include <hydroc/core/system_state.h>
#include <hydroc/radiation/radiation_types.h>
#include "../radiation/radiation_ss_model.h"
#include "../radiation/radiation_ss_fitter.h"
#include <Eigen/Dense>
#include <vector>
#include <memory>

class HydroData;

namespace hydrochrono::hydro {

/**
 * @brief Entry for a single fitted DOF-pair model.
 * 
 * Stores the indices (i,j) and the state-space model for K_ij(t).
 */
struct DofPairModel {
    int i;      ///< Output DOF index (0 to 6*num_bodies - 1)
    int j;      ///< Input DOF index (0 to 6*num_bodies - 1)
    double r2;  ///< Fit quality (R²)
    RadiationStateSpaceModel model;  ///< The fitted model
};

/**
 * @brief Radiation damping force component (state-space approximation).
 * 
 * Computes radiation damping forces using state-space models fitted to
 * RIRF kernels. This provides O(1) per-timestep computation vs O(N) for
 * direct convolution.
 */
class RadiationStateSpaceComponent : public IHydroForceComponent {
public:
    /**
     * @brief Constructor.
     * 
     * Reads RIRF kernels from HydroData, fits state-space models to each
     * DOF-pair kernel, and prepares for runtime computation.
     * 
     * @param file_info Reference to HydroData containing RIRF kernels
     * @param num_bodies Number of bodies in the system
     * @param options State-space fitting options (max_order, r2_threshold)
     */
    RadiationStateSpaceComponent(const HydroData& file_info,
                                  int num_bodies,
                                  const StateSpaceOptions& options);

    /**
     * @brief Get the component type.
     * @return HydroComponentType::Radiation
     */
    HydroComponentType Type() const override { return HydroComponentType::Radiation; }

    /**
     * @brief Compute radiation damping force contribution.
     * 
     * Steps all state-space models forward with current velocities,
     * accumulates forces, and adds to inout_forces per body.
     * 
     * @param state Current system state (velocities needed)
     * @param time Current simulation time
     * @param inout_forces Force vector to add contribution to
     */
    void Compute(const SystemState& state,
                 double time,
                 BodyForces& inout_forces) override;

    /**
     * @brief Reset all internal states to zero.
     * 
     * Call this to restart the simulation without recreating the component.
     */
    void Reset();

    // Diagnostics accessors
    int num_dof_pairs() const { return static_cast<int>(dof_pair_models_.size()); }
    int total_modes() const;
    double min_r2() const;
    double max_r2() const;

private:
    int num_bodies_;
    int num_dofs_;  ///< Total DOFs = 6 * num_bodies
    StateSpaceOptions options_;
    
    std::vector<DofPairModel> dof_pair_models_;  ///< Fitted models for all DOF pairs
    
    double prev_time_;
    bool has_prev_time_;
    
    static constexpr int kDofPerBody = 6;
    
    /**
     * @brief Initialize by fitting models to RIRF kernels.
     * 
     * For each (i,j) DOF pair:
     *   1. Extract K_ij(t) from HydroData
     *   2. Fit using RadiationStateSpaceFitter
     *   3. Create RadiationStateSpaceModel
     *   4. Store in dof_pair_models_
     * 
     * @param file_info HydroData containing RIRF kernels
     */
    void InitializeFromRirf(const HydroData& file_info);
    
    /**
     * @brief Extract velocity vector from system state.
     * 
     * @param state Current system state
     * @return Velocity vector [6*num_bodies] (linear, angular for each body)
     */
    Eigen::VectorXd ExtractVelocityVector(const SystemState& state) const;
    
    /**
     * @brief Accumulate forces into BodyForces structure.
     * 
     * @param f_rad Total radiation force vector [6*num_bodies]
     * @param inout_forces Force output structure (per body)
     */
    void AccumulateForces(const Eigen::VectorXd& f_rad,
                          BodyForces& inout_forces) const;
};

}  // namespace hydrochrono::hydro

#endif  // HYDRO_FORCE_COMPONENTS_RADIATION_SS_COMPONENT_H

