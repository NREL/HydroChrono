// =============================================================================
// HydroChrono VSG Water Surface Visualization
// =============================================================================
// Provides both static and animated water surface rendering for ocean simulation.
//
// AnimatedWaterSurface:
//   Manages a deformable water surface mesh for visualizing ocean waves.
//   Each frame, queries the wave model for elevation eta(x,y,t) at each vertex
//   and updates the GPU buffers accordingly.
//
// CreateStaticWaterPlane:
//   Creates a simple flat water plane for use when no wave model is set.
// =============================================================================
#pragma once

#include <memory>
#include <string>

#include <hydroc/waves/wave_base.h>

#include <chrono/geometry/ChTriangleMeshConnected.h>
#include <chrono/physics/ChBody.h>
#include <chrono_vsg/ChVisualSystemVSG.h>

#include <vsg/all.h>

namespace hydroc {
namespace gui {

// Forward declaration (defined in vsg_gui_component.h).
struct ViewerSettings;

// =============================================================================
// AnimatedWaterSurface
// =============================================================================
/// Manages a deformable water surface mesh for visualizing ocean waves.
///
/// How it works:
/// 1. Creates a flat grid of triangles (64×64 vertices by default).
/// 2. Each frame, queries the wave model for elevation eta(x,y,t) at each vertex.
/// 3. Updates the Z-coordinate of each vertex to match the wave height.
/// 4. Recomputes surface normals for correct lighting.
/// 5. Tells the GPU to re-upload the modified vertex data.
///
/// The mesh uses "soup" format where each triangle has its own 3 vertices
/// (not shared with neighbors), which simplifies GPU buffer updates.
class AnimatedWaterSurface {
  public:
    AnimatedWaterSurface();
    ~AnimatedWaterSurface();

    // Non-copyable, non-movable (owns GPU resources).
    AnimatedWaterSurface(const AnimatedWaterSurface&) = delete;
    AnimatedWaterSurface& operator=(const AnimatedWaterSurface&) = delete;

    /// Reset state for reuse with a new visual system.
    void Reset();

    /// Initialize the mesh geometry and add to the visual system scene.
    /// @note Must be called after ChVisualSystemVSG::Initialize() so the scene exists.
    /// @param vis Pointer to the initialized VSG visual system.
    /// @param resolution Grid resolution (vertices per side). Uses default if 0.
    void Initialize(chrono::vsg3d::ChVisualSystemVSG* vis, int resolution = 0);

    /// Reinitialize with a new grid resolution. Safe to call if already initialized.
    /// @param resolution New grid resolution (vertices per side).
    void Reinitialize(int resolution);

    /// Update vertex Z coordinates based on wave elevation at time t.
    /// If wave is null, keeps surface flat (useful for forced water surface mode).
    /// @param wave The wave model to query for elevation (may be null).
    /// @param t Current simulation time in seconds.
    /// @param settings Optional viewer settings (for scale multiplier and throttle).
    void Update(const std::shared_ptr<WaveBase>& wave, double t,
                const ViewerSettings* settings = nullptr);

    /// Set visibility of the water surface.
    /// @param visible True to show, false to hide (skips rendering).
    void SetVisible(bool visible);

    /// Check if the water surface is currently visible.
    bool IsVisible() const;

    /// Check if the water surface has been initialized.
    bool IsInitialized() const;

    /// Check if initialized for the given visual system.
    /// @param vis The visual system to check against.
    bool IsInitializedFor(chrono::vsg3d::ChVisualSystemVSG* vis) const;

    /// Get current grid resolution.
    int GetGridResolution() const { return current_resolution_; }

    /// Get a human-readable status string (for logging).
    std::string GetStatusString() const;

  private:
    /// Internal initialization with specified resolution.
    void InitializeInternal(chrono::vsg3d::ChVisualSystemVSG* vis, int resolution);

    std::shared_ptr<chrono::ChTriangleMeshConnected> mesh_;
    vsg::ref_ptr<vsg::vec3Array> vsg_vertices_;   ///< VSG vertex positions (dynamic, mesh soup)
    vsg::ref_ptr<vsg::vec3Array> vsg_normals_;    ///< VSG vertex normals (dynamic, mesh soup)
    vsg::ref_ptr<vsg::Node> vsg_node_;            ///< Reference to VSG scene node
    vsg::ref_ptr<vsg::Group> scene_;              ///< Cached scene for add/remove
    chrono::vsg3d::ChVisualSystemVSG* bound_vis_ = nullptr;
    size_t num_triangles_ = 0;                    ///< Number of triangles in mesh
    int current_resolution_ = 0;                  ///< Current grid resolution
    bool initialized_ = false;
    bool visible_ = true;                         ///< Visibility state
    mutable int frame_count_ = 0;                 ///< Frame counter for debug output throttling
    double last_update_time_ = -1.0;              ///< Last update time (for throttling)
};

// =============================================================================
// Static Water Plane Factory
// =============================================================================

/// Creates a visual-only water plane (fixed, non-colliding) at z ~ 0.
/// Used when no wave model is set (static still water).
/// @return A shared pointer to the water plane body (ready to add to system).
std::shared_ptr<chrono::ChBody> CreateStaticWaterPlane();

}  // namespace gui
}  // namespace hydroc

