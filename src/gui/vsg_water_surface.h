// HydroChrono VSG Water Surface Visualization
#pragma once

#include <memory>
#include <string>

#include <hydroc/waves/wave_base.h>

#include <chrono/geometry/ChTriangleMeshConnected.h>
#include <chrono/physics/ChBody.h>
#include <chrono_vsg/ChVisualSystemVSG.h>

#include <vsg/all.h>

#include "vsg_radiation_surface.h"

namespace hydroc {
namespace gui {

// Forward declaration (defined in vsg_gui_component.h).
struct ViewerSettings;

/// Deformable water surface mesh for visualizing ocean waves.
/// Updates vertex Z-coordinates each frame based on wave elevation.
class AnimatedWaterSurface {
  public:
    AnimatedWaterSurface();
    ~AnimatedWaterSurface();

    // Non-copyable, non-movable (owns GPU resources).
    AnimatedWaterSurface(const AnimatedWaterSurface&) = delete;
    AnimatedWaterSurface& operator=(const AnimatedWaterSurface&) = delete;

    /// Reset state for reuse with a new visual system.
    void Reset();

    /// Initialize mesh and add to scene (call after ChVisualSystemVSG::Initialize).
    void Initialize(chrono::vsg3d::ChVisualSystemVSG* vis, int resolution = 0,
                    const ViewerSettings* settings = nullptr);

    /// Reinitialize with new grid resolution or extent.
    void Reinitialize(int resolution, const ViewerSettings* settings = nullptr);

    /// Update vertex Z from wave elevation at time t.
    void Update(const std::shared_ptr<WaveBase>& wave, double t,
                const ViewerSettings* settings = nullptr);

    void SetVisible(bool visible);
    bool IsVisible() const;
    bool IsInitialized() const;
    bool IsInitializedFor(chrono::vsg3d::ChVisualSystemVSG* vis) const;
    int GetGridResolution() const { return current_resolution_; }
    std::string GetStatusString() const;

    RadiationSurfaceViz& GetRadiationViz() { return radiation_viz_; }
    const RadiationSurfaceViz& GetRadiationViz() const { return radiation_viz_; }

    void SetWireframeVisible(bool visible);
    bool IsWireframeVisible() const { return wireframe_visible_; }

  private:
    void InitializeInternal(chrono::vsg3d::ChVisualSystemVSG* vis, int resolution,
                            const ViewerSettings* settings = nullptr);
    void InitializeWireframe();
    void UpdateWireframe();

    std::shared_ptr<chrono::ChTriangleMeshConnected> mesh_;
    vsg::ref_ptr<vsg::vec3Array> vsg_vertices_;
    vsg::ref_ptr<vsg::vec3Array> vsg_normals_;
    vsg::ref_ptr<vsg::vec4Array> vsg_colors_;  ///< Per-vertex colors for height shading
    vsg::ref_ptr<vsg::Node> vsg_node_;
    vsg::ref_ptr<vsg::Group> scene_;
    chrono::vsg3d::ChVisualSystemVSG* bound_vis_ = nullptr;
    size_t num_triangles_ = 0;
    int current_resolution_ = 0;
    bool initialized_ = false;
    bool visible_ = true;
    mutable int frame_count_ = 0;
    double last_update_time_ = -1.0;

    RadiationSurfaceViz radiation_viz_;

    vsg::ref_ptr<vsg::vec3Array> wireframe_vertices_;
    vsg::ref_ptr<vsg::Node> wireframe_node_;
    bool wireframe_visible_ = false;
    bool wireframe_initialized_ = false;

    // Adaptive height-shading range (smoothed min/max eta values).
    float adaptive_eta_min_ = 0.0f;
    float adaptive_eta_max_ = 0.0f;
    bool adaptive_range_initialized_ = false;

    // Custom grid extent. Values <= 0 use defaults from vsg_config.h.
    double grid_width_ = 0.0;   ///< Grid extent in X direction [m]
    double grid_length_ = 0.0;  ///< Grid extent in Y direction [m]
    double grid_center_x_ = 0.0;  ///< Grid center X coordinate [m]
    double grid_center_y_ = 0.0;  ///< Grid center Y coordinate [m]
};

/// Creates a static water plane (visual only, no collision).
std::shared_ptr<chrono::ChBody> CreateStaticWaterPlane();

}  // namespace gui
}  // namespace hydroc

