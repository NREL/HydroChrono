/*********************************************************************
 * @file  mooring_component.cpp
 * @brief Implementation of MooringComponent.
 *********************************************************************/

#include "mooring_component.h"
#include "../mooring/moordyn_wrapper.h"

#include <cassert>
#include <cmath>

namespace hydrochrono::hydro {

MooringComponent::MooringComponent(std::unique_ptr<MoorDynWrapper> wrapper)
    : wrapper_(std::move(wrapper)) {
    assert(wrapper_ && "MooringComponent requires a non-null MoorDynWrapper");
}

void MooringComponent::Compute(const SystemState& state,
                               double time,
                               BodyForces& inout_forces) {
    if (!wrapper_ || !wrapper_->IsInitialized()) {
        return;
    }

    // Tolerance for detecting "same time" (HHT Newton re-evaluations).
    constexpr double kTimeTol = 1.0e-12;
    const bool same_time = has_cached_forces_ &&
                           std::abs(time - last_time_) < kTimeTol;

    if (same_time) {
        ApplyCachedForces(inout_forces);
        return;
    }

    double dt = 0.0;
    if (last_time_ < 0.0) {
        dt = 1.0e-6;
    } else {
        dt = time - last_time_;
    }
    last_time_ = time;

    if (dt <= 0.0) {
        return;
    }

    // Step MoorDyn into a temporary buffer so we can cache the pure
    // mooring contribution separately from the accumulated inout_forces.
    cached_forces_.assign(inout_forces.size(), GeneralizedForce::Zero(6));
    wrapper_->Step(state, time, dt, cached_forces_);
    has_cached_forces_ = true;

    ApplyCachedForces(inout_forces);
}

void MooringComponent::ApplyCachedForces(BodyForces& inout_forces) const {
    for (size_t i = 0; i < cached_forces_.size(); ++i) {
        if (cached_forces_[i].size() == 0) continue;
        if (inout_forces[i].size() < cached_forces_[i].size()) {
            inout_forces[i].setZero(cached_forces_[i].size());
        }
        inout_forces[i] += cached_forces_[i];
    }
}

}  // namespace hydrochrono::hydro
