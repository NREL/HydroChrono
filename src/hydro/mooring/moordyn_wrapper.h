/*********************************************************************
 * @file  moordyn_wrapper.h
 * @brief RAII wrapper around MoorDyn's C API for mooring force computation.
 *
 * MAIN TYPES:
 *   - MoorDynWrapper: Owns a MoorDyn instance, maps SystemState to/from
 *     MoorDyn's flat double arrays, and drives the MoorDyn time stepper.
 *
 * ROLE: Encapsulates all MoorDyn-specific logic so that the rest of
 * HydroChrono only interacts through SystemState and BodyForces.
 *********************************************************************/

#ifndef HYDROC_MOORING_MOORDYN_WRAPPER_H
#define HYDROC_MOORING_MOORDYN_WRAPPER_H

#include <hydroc/core/system_state.h>
#include <hydroc/core/hydro_types.h>

#include <string>
#include <vector>

struct __MoorDyn;
typedef struct __MoorDyn* MoorDyn;

namespace hydrochrono::hydro {

class MoorDynWrapper {
public:
    /**
     * @brief Construct wrapper (does NOT initialize MoorDyn yet).
     *
     * @param input_file Path to MoorDyn input file (e.g. "mooring/lines.txt")
     * @param coupled_body_indices 0-based indices into SystemState::bodies
     *        that correspond to MoorDyn coupled bodies, in the order they
     *        appear in the MoorDyn input file.
     */
    MoorDynWrapper(const std::string& input_file,
                   const std::vector<int>& coupled_body_indices);

    ~MoorDynWrapper();

    MoorDynWrapper(const MoorDynWrapper&) = delete;
    MoorDynWrapper& operator=(const MoorDynWrapper&) = delete;

    /**
     * @brief Create the MoorDyn system and compute initial conditions.
     *
     * Must be called once before the first Step(). Uses the initial body
     * positions/velocities from @p initial_state for the coupled bodies.
     */
    void Initialize(const SystemState& initial_state);

    /**
     * @brief Advance MoorDyn by dt and retrieve mooring forces.
     *
     * @param state   Current system state (positions, velocities)
     * @param time    Current simulation time (seconds)
     * @param dt      Time step (seconds); must be > 0
     * @param[out] out_forces  Forces are ADDED to this vector for the
     *             coupled bodies. Caller must pre-size to num_bodies.
     */
    void Step(const SystemState& state, double time, double dt,
              BodyForces& out_forces);

    bool IsInitialized() const { return initialized_; }

private:
    void PackState(const SystemState& state,
                   std::vector<double>& x,
                   std::vector<double>& xd) const;

    std::string input_file_;
    std::vector<int> coupled_body_indices_;

    MoorDyn system_ = nullptr;
    unsigned int n_coupled_dof_ = 0;
    bool initialized_ = false;
    double moordyn_time_ = 0.0;

    std::vector<double> x_;
    std::vector<double> xd_;
    std::vector<double> f_;
};

}  // namespace hydrochrono::hydro

#endif  // HYDROC_MOORING_MOORDYN_WRAPPER_H
