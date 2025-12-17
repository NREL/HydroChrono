// =============================================================================
// HydroChrono VSG Lighting Utilities
// =============================================================================
// Functions for configuring the VSG scene lighting setup.
// Implements a 3-light studio configuration (key + fill + ambient).
// =============================================================================
#pragma once

#include <chrono_vsg/ChVisualSystemVSG.h>

namespace hydroc {
namespace gui {

/// Adds a fill light to the VSG scene for improved shape readability.
/// The fill light provides soft illumination from the opposite side of the key light,
/// reducing harsh shadows and lifting the dark side of objects.
///
/// @note Must be called AFTER ChVisualSystemVSG::Initialize() since the scene
///       is created during initialization.
///
/// @param vis Pointer to the initialized VSG visual system. Does nothing if null.
void AddFillLight(chrono::vsg3d::ChVisualSystemVSG* vis);

}  // namespace gui
}  // namespace hydroc

