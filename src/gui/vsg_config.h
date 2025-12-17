// =============================================================================
// HydroChrono VSG Configuration Constants
// =============================================================================
// All tuning parameters for the VSG visualization backend in one place.
// These are compile-time constants (constexpr) for zero runtime overhead.
//
// Organized by category:
//   - Appearance (materials, colors)
//   - Water surface (static and animated)
//   - Camera
//   - Lighting
//   - Debug instrumentation
// =============================================================================
#pragma once

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace hydroc {
namespace gui {
namespace vsg_config {

// =============================================================================
// Appearance: Painted Metal (Industrial Yellow)
// =============================================================================
// Base color for rigid bodies - industrial equipment aesthetic.
inline constexpr float kPaintedMetalR = 0.92f;
inline constexpr float kPaintedMetalG = 0.72f;
inline constexpr float kPaintedMetalB = 0.12f;

// Moderate, subtle highlight (painted steel: not mirror, not plastic).
inline constexpr float kSpecular = 0.25f;
inline constexpr float kSpecularExponent = 64.0f;

// Metallic/roughness parameters (painted metal: mostly dielectric coating).
inline constexpr float kRoughness = 0.45f;
inline constexpr float kMetallic = 0.05f;

// Per-body color variation range [0.92, 1.08] for visual distinction.
inline constexpr float kColorVariationMin = 0.92f;
inline constexpr float kColorVariationMax = 1.08f;
inline constexpr int kColorVariationCycle = 8;  // Cycle length for deterministic variation

// =============================================================================
// Appearance: Water Surface
// =============================================================================
// Translucent ocean blue — tropical/deep water aesthetic.
inline constexpr float kWaterR = 0.01f;
inline constexpr float kWaterG = 0.20f;
inline constexpr float kWaterB = 0.35f;

// Translucent surface (allows seeing through to objects below).
inline constexpr float kWaterOpacity = 0.55f;  // 0=fully transparent, 1=fully opaque

// Specular reflection intensity (white highlight on water surface).
inline constexpr float kWaterSpecular = 0.6f;

// Glossy surface: low roughness = more reflective (like calm water).
inline constexpr float kWaterRoughness = 0.05f;
inline constexpr float kWaterMetallic = 0.0f;  // Water is non-metallic (dielectric)

// =============================================================================
// Water Plane Geometry (Static)
// =============================================================================
inline constexpr double kWaterPlaneSize = 100.0;       // meters (large enough to fill most views)
inline constexpr double kWaterPlaneThickness = 0.02;   // meters (thin slab)
inline constexpr double kWaterPlaneZ = -0.05;          // below z=0 to avoid z-fighting

// =============================================================================
// Water Grid (Animated Surface)
// =============================================================================
inline constexpr double kWaterGridSize = 100.0;        // meters (grid extent in x and y)
inline constexpr int kWaterGridResolution = 100;        // vertices per side (64×64)
inline constexpr double kWaterGridSpacing = kWaterGridSize / (kWaterGridResolution - 1);

// =============================================================================
// Camera
// =============================================================================
// Deterministic side-on view: eye on negative Y axis, slightly above, looking at origin.
inline constexpr double kCameraDistance = 60.0;        // meters along -Y axis
inline constexpr double kCameraHeight = 40.0;          // meters above z=0

// =============================================================================
// Lighting: 3-Light Studio Setup (Key + Fill + Ambient)
// =============================================================================
// Key light (main sun) — primary illumination from upper-front.
inline constexpr float kKeyIntensity = 1.0f;
inline constexpr double kKeyAzimuth = 0.4 * M_PI;      // ~72° from front
inline constexpr double kKeyElevation = M_PI / 5.0;    // ~36° above horizon

// Fill light (soft opposite) — reduces harsh shadows, lifts dark side.
inline constexpr float kFillIntensity = 0.35f;
inline constexpr double kFillAzimuth = kKeyAzimuth + M_PI;  // opposite side
inline constexpr double kFillElevation = M_PI / 8.0;        // ~22.5° above horizon

// Note: Chrono-VSG includes a built-in ambient at 10% of key intensity.
// Shadows remain OFF (water plane lacks per-object shadow control).

// =============================================================================
// Wireframe Overlay
// =============================================================================
// Enable wireframe rendering for improved shape readability.
// Default OFF: solid shaded rendering displays materials/colors properly.
inline constexpr bool kEnableWireframe = false;

// =============================================================================
// Debug Instrumentation
// =============================================================================
// Set to true to diagnose flat free-surface issues.
inline constexpr bool kDebugWaveSurface = false;       // OFF by default
inline constexpr int kDebugPrintEveryNFrames = 60;     // Print every N frames to reduce spam
inline constexpr bool kDebugSpikeVertex = false;       // OFF by default; verify VSG pipeline
inline constexpr bool kForceWaterSurface = false;      // OFF by default; force water even without waves

}  // namespace vsg_config
}  // namespace gui
}  // namespace hydroc

