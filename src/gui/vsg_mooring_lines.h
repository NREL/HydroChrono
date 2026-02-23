/// @file vsg_mooring_lines.h
/// @brief VSG mooring-line tube-mesh renderer for HydroChrono.
#pragma once

#include <hydroc/gui/mooring_viz_data.h>

#include <chrono_vsg/ChVisualSystemVSG.h>

#include <vsg/all.h>

#include <vector>

namespace hydroc {
namespace gui {

/// Renders MoorDyn mooring lines as cylindrical tube meshes in the VSG scene.
///
/// Each line is tessellated once during Initialize() and then cheaply
/// repositioned every frame via direct writes into the VSG vertex buffer,
/// followed by a `dirty()` call that triggers GPU re-upload.
class MooringLinesViz {
  public:
    MooringLinesViz();
    ~MooringLinesViz();

    MooringLinesViz(const MooringLinesViz&) = delete;
    MooringLinesViz& operator=(const MooringLinesViz&) = delete;

    /// Build tube meshes from initial line data and add to the VSG scene.
    /// Call once after ChVisualSystemVSG::Initialize().
    ///
    /// @param vis            Pointer to the active visual system (must outlive
    ///                       this object).
    /// @param initial_lines  One entry per mooring line with initial node
    ///                       positions and endpoint types.
    /// @param parent_group   Optional scene-graph group.  When provided,
    ///                       mooring geometry is added here instead of the root
    ///                       VSG scene.  Use this to guarantee render order
    ///                       (e.g. render mooring lines before the transparent
    ///                       water surface).
    void Initialize(chrono::vsg3d::ChVisualSystemVSG* vis,
                    const std::vector<MooringLineVizData>& initial_lines,
                    vsg::ref_ptr<vsg::Group> parent_group = {});

    /// Reposition tube vertices to match current node positions.
    void Update(const std::vector<MooringLineVizData>& lines);

    bool IsInitialized() const { return initialized_; }
    bool IsInitializedFor(chrono::vsg3d::ChVisualSystemVSG* vis) const {
        return initialized_ && bound_vis_ == vis;
    }

  private:
    static constexpr int kSides = 8;
    static constexpr double kTubeRadius = 0.15;
    static constexpr double kEndpointRadius = 0.5;
    static constexpr double kNodeMarkerRadius = 0.20;

    /// A PBR sphere used to mark an endpoint or intermediate node.
    struct PointMarker {
        vsg::ref_ptr<vsg::MatrixTransform> transform;
        vsg::ref_ptr<vsg::Node> node;
    };

    /// Per-line rendering state: VSG node, extracted vertex/normal buffers
    /// (for dynamic updates), and small sphere markers at key positions.
    struct LineMesh {
        vsg::ref_ptr<vsg::vec3Array> vertices;
        vsg::ref_ptr<vsg::vec3Array> normals;
        vsg::ref_ptr<vsg::Node> node;
        size_t num_nodes = 0;
        PointMarker start_marker;
        PointMarker end_marker;
        std::vector<PointMarker> node_markers;
    };

    void BuildTubeMesh(const MooringLineVizData& line_data, LineMesh& lm,
                       chrono::vsg3d::ChVisualSystemVSG* vis);
    void UpdateTubeMesh(const MooringLineVizData& line_data, LineMesh& lm);

    PointMarker CreateSphereMarker(chrono::vsg3d::ChVisualSystemVSG* vis,
                                   double x, double y, double z,
                                   double radius, int point_type);

    static void UpdateMarkerPosition(PointMarker& marker,
                                     double x, double y, double z,
                                     double radius);

    std::vector<LineMesh> line_meshes_;
    vsg::ref_ptr<vsg::Group> scene_;
    chrono::vsg3d::ChVisualSystemVSG* bound_vis_ = nullptr;
    bool initialized_ = false;
};

}  // namespace gui
}  // namespace hydroc
