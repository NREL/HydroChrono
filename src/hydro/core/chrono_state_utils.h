/*********************************************************************
 * @file  chrono_state_utils.h
 * @brief Utilities for converting between Chrono and SystemState representations.
 *********************************************************************/

#ifndef HYDRO_CORE_CHRONO_STATE_UTILS_H
#define HYDRO_CORE_CHRONO_STATE_UTILS_H

#include <hydroc/system_state.h>
#include <chrono/physics/ChBody.h>
#include <memory>
#include <vector>

namespace hydrochrono::hydro {

/**
 * @brief Build SystemState from Chrono bodies.
 * 
 * Extracts pose and velocity data from Chrono bodies into the SystemState structure.
 * 
 * @param bodies Vector of Chrono body pointers
 * @param out_state Output SystemState (will be populated)
 */
void BuildSystemStateFromChronoBodies(
    const std::vector<std::shared_ptr<chrono::ChBody>>& bodies,
    SystemState& out_state);

}  // namespace hydrochrono::hydro

#endif  // HYDRO_CORE_CHRONO_STATE_UTILS_H

