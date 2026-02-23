/*********************************************************************
 * @file  mooring_component.h
 * @brief Mooring force component backed by MoorDyn.
 *********************************************************************/

#ifndef HYDRO_FORCE_COMPONENTS_MOORING_COMPONENT_H
#define HYDRO_FORCE_COMPONENTS_MOORING_COMPONENT_H

#include <hydroc/core/force_component.h>
#include <hydroc/core/system_state.h>
#include <hydroc/gui/mooring_viz_data.h>
#include <memory>
#include <vector>

namespace hydrochrono::hydro {

class MoorDynWrapper;

/**
 * @brief IHydroForceComponent implementation that delegates to MoorDyn
 *        via MoorDynWrapper for mooring line force computation.
 *
 * Caches the most recent MoorDyn forces so that repeated evaluations
 * at the same time (e.g. HHT Newton iterations) return consistent
 * values without re-stepping MoorDyn.
 */
class MooringComponent : public IHydroForceComponent {
public:
    explicit MooringComponent(std::unique_ptr<MoorDynWrapper> wrapper);

    HydroComponentType Type() const override {
        return HydroComponentType::Mooring;
    }

    void Compute(const SystemState& state,
                 double time,
                 BodyForces& inout_forces) override;

    /// Return current node positions for all mooring lines (forwarded from MoorDynWrapper).
    std::vector<hydroc::gui::MooringLineVizData> GetMooringLineStates() const;

private:
    void ApplyCachedForces(BodyForces& inout_forces) const;

    std::unique_ptr<MoorDynWrapper> wrapper_;
    double last_time_ = -1.0;
    BodyForces cached_forces_;
    bool has_cached_forces_ = false;
};

}  // namespace hydrochrono::hydro

#endif  // HYDRO_FORCE_COMPONENTS_MOORING_COMPONENT_H
