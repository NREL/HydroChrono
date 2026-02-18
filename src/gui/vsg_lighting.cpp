// =============================================================================
// HydroChrono VSG Lighting Utilities - Implementation
// =============================================================================
#include "vsg_lighting.h"

#include "vsg_config.h"

#include <cmath>

#include <vsg/all.h>

namespace hydroc {
namespace gui {

using namespace vsg_config;

void AddFillLight(chrono::vsg3d::ChVisualSystemVSG* vis) {
    if (!vis) {
        return;
    }

    auto scene = vis->GetVSGScene();
    if (!scene) {
        return;
    }

    // Compute fill light direction (Z-up coordinate system).
    const double se = std::sin(kFillElevation);
    const double ce = std::cos(kFillElevation);
    const double sa = std::sin(kFillAzimuth);
    const double ca = std::cos(kFillAzimuth);

    auto fill_light = vsg::DirectionalLight::create();
    fill_light->name = "fill light";
    fill_light->color.set(1.0f, 1.0f, 1.0f);
    fill_light->intensity = kFillIntensity;
    fill_light->direction.set(
        static_cast<float>(-ce * ca),
        static_cast<float>(-ce * sa),
        static_cast<float>(-se));

    scene->addChild(fill_light);
}

}  // namespace gui
}  // namespace hydroc

