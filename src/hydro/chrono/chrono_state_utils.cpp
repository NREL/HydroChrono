/*********************************************************************
 * @file  chrono_state_utils.cpp
 * @brief Implementation of Chrono-SystemState conversion utilities.
 *
 * Part of the Chrono coupling layer (src/hydro/chrono/).
 *********************************************************************/

#include "chrono_state_utils.h"

namespace hydrochrono::hydro::chrono_coupling {

void BuildSystemStateFromChronoBodies(
    const std::vector<std::shared_ptr<chrono::ChBody>>& bodies,
    SystemState& out_state) {
    out_state.bodies.clear();
    out_state.bodies.reserve(bodies.size());
    
    for (const auto& body : bodies) {
        BodyState body_state;
        
        // Extract position and orientation
        const chrono::ChVector3d position_world = body->GetPos();
        const chrono::ChVector3d rotation_rpy    = body->GetRot().GetCardanAnglesXYZ();
        
        body_state.position[0] = position_world.x();
        body_state.position[1] = position_world.y();
        body_state.position[2] = position_world.z();
        
        body_state.orientation_rpy[0] = rotation_rpy.x();
        body_state.orientation_rpy[1] = rotation_rpy.y();
        body_state.orientation_rpy[2] = rotation_rpy.z();
        
        // Extract velocities
        const auto linear_velocity_world  = body->GetPosDt();
        const auto angular_velocity_world = body->GetAngVelParent();
        
        body_state.linear_velocity[0] = linear_velocity_world.x();
        body_state.linear_velocity[1] = linear_velocity_world.y();
        body_state.linear_velocity[2] = linear_velocity_world.z();
        
        body_state.angular_velocity[0] = angular_velocity_world.x();
        body_state.angular_velocity[1] = angular_velocity_world.y();
        body_state.angular_velocity[2] = angular_velocity_world.z();
        
        out_state.bodies.push_back(body_state);
    }
}

}  // namespace hydrochrono::hydro::chrono_coupling

