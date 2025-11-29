/*********************************************************************
 * @file  chrono_state_utils.h
 * @brief Utilities for converting between Chrono and SystemState representations.
 *
 * This file is part of the Chrono coupling layer (src/hydro/chrono/).
 * It provides the bridge between Chrono physics types (ChBody) and the
 * Chrono-free core hydrodynamics types (SystemState).
 *
 * ARCHITECTURE NOTE:
 * - Core hydrodynamics (src/hydro/core/) remains Chrono-free.
 * - This utility lives in the coupling layer because it depends on Chrono types.
 * - The core SystemState types are defined in include/hydroc/core/system_state.h.
 *********************************************************************/

#ifndef HYDRO_CHRONO_CHRONO_STATE_UTILS_H
#define HYDRO_CHRONO_CHRONO_STATE_UTILS_H

#include <hydroc/core/system_state.h>
#include <chrono/physics/ChBody.h>
#include <memory>
#include <vector>

namespace hydrochrono::hydro::chrono_coupling {

/**
 * @brief Build SystemState from Chrono bodies.
 * 
 * Extracts pose and velocity data from Chrono bodies into the SystemState structure.
 * This function bridges the Chrono physics layer with the Chrono-free core.
 * 
 * @param bodies Vector of Chrono body pointers
 * @param out_state Output SystemState (will be populated)
 */
void BuildSystemStateFromChronoBodies(
    const std::vector<std::shared_ptr<chrono::ChBody>>& bodies,
    SystemState& out_state);

}  // namespace hydrochrono::hydro::chrono_coupling

#endif  // HYDRO_CHRONO_CHRONO_STATE_UTILS_H

