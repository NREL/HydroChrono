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

