// =============================================================================
// HydroChrono VSG Material Utilities
// =============================================================================
// Factory functions for creating and applying visual materials to bodies.
// Provides consistent industrial appearance across the visualization.
// =============================================================================
#pragma once

#include <chrono/assets/ChVisualMaterial.h>
#include <chrono/physics/ChBody.h>

namespace hydroc {
namespace gui {

/// Creates the base painted-metal material (industrial yellow).
/// @return A ChVisualMaterial configured with painted metal appearance.
chrono::ChVisualMaterial MakePaintedMetalBase();

/// Creates a painted-metal variant with subtle brightness variation for body distinction.
/// The variation is deterministic based on the body index, cycling through a fixed range.
/// @param index Body index used for deterministic color variation.
/// @return A ChVisualMaterial with slight color variation from the base.
chrono::ChVisualMaterial MakePaintedMetalVariant(int index);

/// Creates the water surface material (translucent ocean blue).
/// @return A ChVisualMaterial configured for water appearance with transparency.
chrono::ChVisualMaterial MakeWaterSurfaceMaterial();

/// Applies a material to all visual shapes attached to a body.
/// Replaces existing materials or adds the material if none exist.
/// @param body The body whose visual shapes will receive the material.
/// @param mat The material to apply.
void ApplyMaterialToAllVisualShapes(chrono::ChBody& body, const chrono::ChVisualMaterial& mat);

}  // namespace gui
}  // namespace hydroc

