// =============================================================================
// HydroChrono VSG GUI Component
// =============================================================================
// ImGui overlay component for the VSG visualization backend.
// Provides the play/pause button and other HydroChrono-specific UI elements.
// =============================================================================
#pragma once

#include <chrono_vsg/ChGuiComponentVSG.h>
#include <chrono_vsg/ChVisualSystemVSG.h>

namespace hydroc {
namespace gui {

// =============================================================================
// ViewerSettings: Runtime-adjustable water visualization settings
// =============================================================================
// Allows tuning water appearance without recompiling.
// Owned by GUIImplVSG, passed by pointer to GUI component and water surface.
struct ViewerSettings {
    bool show_water = true;           ///< Toggle water surface visibility
    float wave_visual_scale = 1.0f;   ///< Multiplier for wave elevation display
    int grid_resolution = 64;         ///< Vertices per side (32/64/96/128)
    int update_hz = 30;               ///< Max update rate in Hz (10-60)
    bool show_water_status = false;   ///< Show status line in overlay

    // Allowed discrete values for grid resolution.
    static constexpr int kResolutionOptions[] = {32, 64, 96, 128};
    static constexpr int kResolutionCount = 4;

    // Track if resolution changed (requires mesh rebuild).
    bool resolution_changed = false;

    // =========================================================================
    // Radiation Visualization (Tier 0 - Approximate, Visual Only)
    // =========================================================================
    // Simple visualization of radiated waves from moving bodies.
    // NOT physically accurate - purely for visual feedback.
    //
    // Most parameters are auto-derived from physics:
    //   - Wavelength derived from wave speed using deep water dispersion
    //   - Amplitude scales with body velocity and size
    //   - Only "Visual Scale" is user-adjustable (for inspection)
    bool show_radiation_viz = false;       ///< Enable radiated wave visualization
    float radiation_visual_scale = 1.0f;   ///< Visual amplification (range 0.1-5x, for inspection only)

    // Water surface grid overlay for better visibility.
    bool show_water_grid = false;          ///< Show wireframe grid on water surface
};

/// ImGui component for HydroChrono visualization overlay.
/// Displays a play/pause button and viewer settings panel.
class HydroChronoGuiComponent : public chrono::vsg3d::ChGuiComponentVSG {
  public:
    /// Construct the GUI component.
    /// @param vsys Pointer to the VSG visual system (for potential future use).
    /// @param button_pressed Reference to the simulation running state (toggled by button).
    /// @param settings Pointer to viewer settings (owned by GUIImplVSG).
    HydroChronoGuiComponent(chrono::vsg3d::ChVisualSystemVSG* vsys, bool& button_pressed,
                            ViewerSettings* settings = nullptr);

    /// Render the ImGui overlay.
    /// Called each frame by the VSG rendering loop.
    void render(vsg::CommandBuffer& cb) override;

  private:
    chrono::vsg3d::ChVisualSystemVSG* vsys_;
    bool& pressed_;
    ViewerSettings* settings_;  ///< Viewer settings (may be null)
};

}  // namespace gui
}  // namespace hydroc

