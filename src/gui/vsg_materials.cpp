// =============================================================================
// HydroChrono VSG Material Utilities - Implementation
// =============================================================================
#include "vsg_materials.h"

#include "vsg_config.h"

#include <algorithm>

#include <chrono/core/ChTypes.h>

namespace hydroc {
namespace gui {

using namespace vsg_config;

chrono::ChVisualMaterial MakePaintedMetalBase() {
    chrono::ChVisualMaterial mat;

    mat.SetDiffuseColor(chrono::ChColor(kPaintedMetalR, kPaintedMetalG, kPaintedMetalB));
    mat.SetSpecularColor(chrono::ChColor(kSpecular, kSpecular, kSpecular));
    mat.SetSpecularExponent(kSpecularExponent);
    mat.SetRoughness(kRoughness);
    mat.SetMetallic(kMetallic);
    mat.SetOpacity(1.0f);

    return mat;
}

chrono::ChVisualMaterial MakePaintedMetalVariant(int index) {
    // Compute variation factor in [kColorVariationMin, kColorVariationMax].
    // Cycle through values to keep variation subtle and deterministic.
    const float t = static_cast<float>(index % kColorVariationCycle) /
                    static_cast<float>(kColorVariationCycle - 1);
    const float factor = kColorVariationMin + t * (kColorVariationMax - kColorVariationMin);

    // Apply factor to base color, clamping to valid range.
    const float r = std::min(1.0f, kPaintedMetalR * factor);
    const float g = std::min(1.0f, kPaintedMetalG * factor);
    const float b = std::min(1.0f, kPaintedMetalB * factor);

    chrono::ChVisualMaterial mat;
    mat.SetDiffuseColor(chrono::ChColor(r, g, b));
    mat.SetSpecularColor(chrono::ChColor(kSpecular, kSpecular, kSpecular));
    mat.SetSpecularExponent(kSpecularExponent);
    mat.SetRoughness(kRoughness);
    mat.SetMetallic(kMetallic);
    mat.SetOpacity(1.0f);

    return mat;
}

chrono::ChVisualMaterial MakeWaterSurfaceMaterial() {
    chrono::ChVisualMaterial mat;

    mat.SetDiffuseColor(chrono::ChColor(kWaterR, kWaterG, kWaterB));
    mat.SetSpecularColor(chrono::ChColor(kWaterSpecular, kWaterSpecular, kWaterSpecular));
    mat.SetRoughness(kWaterRoughness);
    mat.SetMetallic(kWaterMetallic);
    mat.SetOpacity(kWaterOpacity);

    return mat;
}

void ApplyMaterialToAllVisualShapes(chrono::ChBody& body, const chrono::ChVisualMaterial& mat) {
    auto model = body.GetVisualModel();
    if (!model) {
        return;
    }

    auto mat_ptr = chrono_types::make_shared<chrono::ChVisualMaterial>(mat);

    const unsigned int num_shapes = model->GetNumShapes();
    for (unsigned int i = 0; i < num_shapes; ++i) {
        auto shape = model->GetShape(i);
        if (!shape) {
            continue;
        }

        const unsigned int num_materials = shape->GetNumMaterials();
        if (num_materials == 0) {
            shape->AddMaterial(mat_ptr);
            continue;
        }

        for (unsigned int m = 0; m < num_materials; ++m) {
            shape->SetMaterial(m, mat_ptr);
        }
    }
}

}  // namespace gui
}  // namespace hydroc

