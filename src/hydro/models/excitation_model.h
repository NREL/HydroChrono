/*********************************************************************
 * @file  excitation_model.h
 * @brief Wave excitation force model.
 *********************************************************************/

#ifndef HYDRO_MODELS_EXCITATION_MODEL_H
#define HYDRO_MODELS_EXCITATION_MODEL_H

#include "../core/hydro_force_model.h"
#include "../core/system_state.h"
#include "../waves/wave_base.h"
#include <memory>

namespace hydrochrono::hydro {

/**
 * @brief Wave excitation force model.
 * 
 * Computes wave excitation forces by calling WaveBase::GetForceAtTime()
 * and converting the flat 6N vector to per-body forces.
 */
class ExcitationModel : public IHydroForceModel {
public:
    /**
     * @brief Constructor.
     * 
     * @param waves Shared pointer to WaveBase instance for force computation
     * @param num_bodies Number of bodies in the system
     */
    ExcitationModel(std::shared_ptr<WaveBase> waves, int num_bodies);

    /**
     * @brief Compute wave excitation force contribution.
     * 
     * Calls wave->GetForceAtTime(time) to get flat 6N force vector,
     * converts to per-body forces, and adds to inout_forces.
     * 
     * @param state Current system state (not used for excitation, but required by interface)
     * @param time Current simulation time
     * @param inout_forces Force vector to add contribution to (one GeneralizedForce per body)
     */
    void Compute(const SystemState& state,
                double time,
                BodyForces& inout_forces) override;

private:
    std::shared_ptr<WaveBase> waves_;
    int num_bodies_;
    
    static constexpr int kDofPerBody = 6;
};

}  // namespace hydrochrono::hydro

#endif  // HYDRO_MODELS_EXCITATION_MODEL_H

