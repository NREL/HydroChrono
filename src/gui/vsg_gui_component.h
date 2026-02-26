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

class MooringLinesViz;  // forward declaration

/// Runtime-adjustable water visualization settings.
///
/// Allows tuning water appearance without recompiling. Owned by GUIImplVSG
/// and passed by pointer to GUI components and the water surface renderer.
struct ViewerSettings {
    // --- Display Options ---
    bool show_water = true;           ///< Toggle water surface visibility
    bool show_water_grid = false;     ///< Show wireframe overlay on water surface
    bool show_water_status = false;   ///< Show status line in overlay

    // --- Wave Rendering ---
    float wave_visual_scale = 1.0f;   ///< Multiplier for wave elevation display
    int update_hz = 30;               ///< Maximum update rate [Hz], range 10-60

    // --- Grid Resolution ---
    int grid_resolution = 64;         ///< Vertices per side (32, 64, 96, or 128)
    bool resolution_changed = false;  ///< Triggers mesh rebuild when resolution changes

    static constexpr int kResolutionOptions[] = {32, 64, 96, 128};
    static constexpr int kResolutionCount = 4;

    // --- Grid Extent ---
    // Values <= 0 use defaults (100m x 100m centered at origin).
    double grid_width = 0.0;          ///< Grid extent in X direction [m]
    double grid_length = 0.0;         ///< Grid extent in Y direction [m]
    double grid_center_x = 0.0;       ///< Grid center X coordinate [m]
    double grid_center_y = 0.0;       ///< Grid center Y coordinate [m]
    bool grid_extent_changed = false; ///< Triggers mesh rebuild when extent changes

    // --- Radiation Visualization (Approximate) ---
    // Visual-only representation of radiated waves from moving bodies.
    // NOT physically accurate; wavelength derived from dispersion relation,
    // amplitude scales with body velocity and size.
    bool show_radiation_viz = false;       ///< Enable radiated wave visualization
    float radiation_visual_scale = 1.0f;   ///< Visual amplification [0.1x - 5x]

    // --- Mooring Colour Mapping ---
    bool show_mooring_colors = true;       ///< Colour lines by tension magnitude
    bool mooring_range_locked = false;     ///< Freeze the adaptive min/max range
};

/// ImGui component for HydroChrono visualization overlay.
/// Displays a play/pause button and viewer settings panel.
class HydroChronoGuiComponent : public chrono::vsg3d::ChGuiComponentVSG {
  public:
    /// Construct the GUI component.
    /// @param vsys Pointer to the VSG visual system (for potential future use).
    /// @param button_pressed Reference to the simulation running state (toggled by button).
    /// @param settings Pointer to viewer settings (owned by GUIImplVSG).
    /// @param mooring_viz Pointer to the mooring renderer (for reading the
    ///                    adaptive scalar range for the colour bar).
    HydroChronoGuiComponent(chrono::vsg3d::ChVisualSystemVSG* vsys, bool& button_pressed,
                            ViewerSettings* settings = nullptr,
                            MooringLinesViz* mooring_viz = nullptr);

    /// Render the ImGui overlay.
    /// Called each frame by the VSG rendering loop.
    void render(vsg::CommandBuffer& cb) override;

  private:
    void RenderMooringPanel();
    void RenderColorBar();

    chrono::vsg3d::ChVisualSystemVSG* vsys_;
    bool& pressed_;
    ViewerSettings* settings_;
    MooringLinesViz* mooring_viz_;
};

}  // namespace gui
}  // namespace hydroc

