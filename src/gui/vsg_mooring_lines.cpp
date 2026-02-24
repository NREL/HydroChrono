/// @file vsg_mooring_lines.cpp
/// @brief Implementation of MooringLinesViz tube-mesh renderer.
#include "vsg_mooring_lines.h"

#include "vsg_config.h"

#include <chrono/assets/ChVisualMaterial.h>
#include <chrono/core/ChTypes.h>
#include <chrono/geometry/ChTriangleMeshConnected.h>
#include <chrono_vsg/utils/ChShapeBuilderVSG.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace hydroc {
namespace gui {

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// VSG visitor that locates the first `vsg::vec3Array` inside a
/// `vsg::VertexIndexDraw` node at array slot @p N (0 = positions, 1 = normals).
template <int N>
class FindVec3Buffer : public vsg::Visitor {
  public:
    void apply(vsg::Object& object) override { object.traverse(*this); }
    void apply(vsg::VertexIndexDraw& vid) override {
        if (static_cast<int>(vid.arrays.size()) <= N) return;
        vid.arrays[N]->data->accept(*this);
    }
    void apply(vsg::vec3Array& arr) override {
        if (!result_) result_ = &arr;
    }
    vsg::ref_ptr<vsg::vec3Array> Result() { return result_; }

  private:
    vsg::ref_ptr<vsg::vec3Array> result_;
};

/// Same as FindVec3Buffer but for `vsg::vec4Array` (vertex colours at slot 3).
template <int N>
class FindVec4Buffer : public vsg::Visitor {
  public:
    void apply(vsg::Object& object) override { object.traverse(*this); }
    void apply(vsg::VertexIndexDraw& vid) override {
        if (static_cast<int>(vid.arrays.size()) <= N) return;
        vid.arrays[N]->data->accept(*this);
    }
    void apply(vsg::vec4Array& arr) override {
        if (!result_) result_ = &arr;
    }
    vsg::ref_ptr<vsg::vec4Array> Result() { return result_; }

  private:
    vsg::ref_ptr<vsg::vec4Array> result_;
};

/// Finite-difference tangent (central where possible, one-sided at ends).
static vsg::dvec3 Tangent(const std::vector<std::array<double, 3>>& pts,
                          size_t i) {
    if (pts.size() < 2) return {0.0, 0.0, 1.0};

    size_t a = 0, b = 0;
    if (i == 0) {
        a = 0;
        b = 1;
    } else if (i >= pts.size() - 1) {
        a = pts.size() - 2;
        b = pts.size() - 1;
    } else {
        a = i - 1;
        b = i + 1;
    }

    vsg::dvec3 t(pts[b][0] - pts[a][0],
                 pts[b][1] - pts[a][1],
                 pts[b][2] - pts[a][2]);
    double len = std::sqrt(t.x * t.x + t.y * t.y + t.z * t.z);
    if (len < 1e-12) return {0.0, 0.0, 1.0};
    return {t.x / len, t.y / len, t.z / len};
}

/// Build an orthonormal frame (u, v) perpendicular to @p tang.
static void PerpFrame(const vsg::dvec3& tang, vsg::dvec3& u, vsg::dvec3& v) {
    const vsg::dvec3 ref =
        (std::abs(tang.z) < 0.9) ? vsg::dvec3(0, 0, 1) : vsg::dvec3(1, 0, 0);

    // u = ref x tang  (cross product)
    u = {ref.y * tang.z - ref.z * tang.y,
         ref.z * tang.x - ref.x * tang.z,
         ref.x * tang.y - ref.y * tang.x};
    double len = std::sqrt(u.x * u.x + u.y * u.y + u.z * u.z);
    if (len < 1e-12) {
        u = {1.0, 0.0, 0.0};
        len = 1.0;
    }
    u = {u.x / len, u.y / len, u.z / len};

    // v = tang x u
    v = {tang.y * u.z - tang.z * u.y,
         tang.z * u.x - tang.x * u.z,
         tang.x * u.y - tang.y * u.x};
}

/// Map a MoorDyn point type (or kIntermediateNode) to a display colour.
static chrono::ChColor PointTypeColor(int point_type) {
    using PT = MooringPointType;
    switch (static_cast<PT>(point_type)) {
        case PT::kFixed:   return {0.6f, 0.15f, 0.1f};   // anchor -- dark rust
        case PT::kCoupled: return {0.2f, 0.4f, 0.7f};    // fairlead -- steel blue
        case PT::kFree:    return {0.9f, 0.6f, 0.1f};    // buoy/free -- orange
    }
    return {0.05f, 0.05f, 0.05f};  // intermediate node -- near black
}

/// Attempt to match the Google Turbo colour map (Mikhailov, 2019).
/// Input @p t is clamped to [0, 1].  Returns an opaque RGBA colour.
static vsg::vec4 TurboColormap(float t) {
    t = std::clamp(t, 0.0f, 1.0f);

    // Polynomial coefficients fitted to the 256-entry Turbo LUT.
    const float r = std::clamp(
        0.13572138f + t * (4.61539260f + t * (-42.66032258f +
        t * (132.13108234f + t * (-152.54895899f + t * 59.28637943f)))),
        0.0f, 1.0f);
    const float g = std::clamp(
        0.09140261f + t * (2.26418794f + t * (4.11868525f +
        t * (-44.58319668f + t * (70.41698018f - t * 33.26974748f)))),
        0.0f, 1.0f);
    const float b = std::clamp(
        0.10667330f + t * (12.75191895f + t * (-60.25290628f +
        t * (109.04872043f + t * (-89.38040853f + t * 27.13073700f)))),
        0.0f, 1.0f);

    return {r, g, b, 1.0f};
}

static const vsg::vec4 kDefaultCableColorVec4{0.55f, 0.45f, 0.35f, 1.0f};

// ---------------------------------------------------------------------------
// Lifetime
// ---------------------------------------------------------------------------

MooringLinesViz::MooringLinesViz() = default;
MooringLinesViz::~MooringLinesViz() = default;

void MooringLinesViz::Initialize(
    chrono::vsg3d::ChVisualSystemVSG* vis,
    const std::vector<MooringLineVizData>& initial_lines,
    vsg::ref_ptr<vsg::Group> parent_group) {
    if (!vis || initial_lines.empty()) return;
    if (initialized_ && bound_vis_ == vis) return;

    bound_vis_ = vis;
    scene_ = parent_group ? parent_group : vis->GetVSGScene();
    if (!scene_) return;

    line_meshes_.resize(initial_lines.size());
    for (size_t li = 0; li < initial_lines.size(); ++li) {
        if (initial_lines[li].node_positions.size() < 2) continue;
        BuildTubeMesh(initial_lines[li], line_meshes_[li], vis);
    }

    initialized_ = true;
    std::cout << "[MooringViz] Initialized " << initial_lines.size()
              << " mooring line(s)\n";
}

void MooringLinesViz::Update(const std::vector<MooringLineVizData>& lines,
                             bool color_enabled, bool range_locked) {
    if (!initialized_) return;
    const size_t count = std::min(lines.size(), line_meshes_.size());

    // Adaptive range: min fixed at 0 (all fields are non-negative magnitudes),
    // max uses high-water-mark behaviour -- snaps to new peaks instantly but
    // decays very slowly so the colour scale stays stable.
    if (color_enabled && !range_locked) {
        float frame_max = 0.0f;
        for (size_t li = 0; li < count; ++li) {
            for (double s : lines[li].node_tensions) {
                frame_max = std::max(frame_max, static_cast<float>(s));
            }
        }

        const float padded_max = frame_max * (1.0f + kRangePadding);

        if (!adaptive_initialized_) {
            adaptive_max_ = std::max(padded_max, 1.0f);
            adaptive_initialized_ = true;
            hold_frames_remaining_ = kRangeHoldFrames;
        } else if (padded_max > adaptive_max_) {
            adaptive_max_ = padded_max;
            hold_frames_remaining_ = kRangeHoldFrames;
        } else if (hold_frames_remaining_ > 0) {
            --hold_frames_remaining_;
        } else {
            adaptive_max_ += kRangeDecayAlpha * (padded_max - adaptive_max_);
        }

        adaptive_min_ = 0.0f;

        if (adaptive_max_ < 1e-6f) {
            adaptive_max_ = 1.0f;
        }
    }

    for (size_t li = 0; li < count; ++li) {
        if (lines[li].node_positions.size() < 2) continue;
        if (!line_meshes_[li].vertices) continue;
        UpdateTubeMesh(lines[li], line_meshes_[li],
                       color_enabled, adaptive_min_, adaptive_max_);
    }
}

// ---------------------------------------------------------------------------

void MooringLinesViz::BuildTubeMesh(const MooringLineVizData& line_data,
                                    LineMesh& lm,
                                    chrono::vsg3d::ChVisualSystemVSG* vis) {
    auto shape_builder = vis->GetVSGShapeBuilder();
    if (!shape_builder) return;

    const auto& pts = line_data.node_positions;
    const size_t num_nodes = pts.size();
    lm.num_nodes = num_nodes;

    // -- Tessellate a closed tube (N rings of kSides verts + 2 cap centres). --

    auto mesh = chrono_types::make_shared<chrono::ChTriangleMeshConnected>();
    auto& verts      = mesh->m_vertices;
    auto& mesh_faces = mesh->m_face_v_indices;

    const size_t ring_verts  = num_nodes * static_cast<size_t>(kSides);
    const size_t total_verts = ring_verts + 2;  // +2 cap centres
    verts.resize(total_verts);

    // Pre-compute trig table for the cross-section.
    std::vector<double> cos_table(kSides);
    std::vector<double> sin_table(kSides);
    for (int s = 0; s < kSides; ++s) {
        const double angle = 2.0 * M_PI * s / kSides;
        cos_table[s] = std::cos(angle);
        sin_table[s] = std::sin(angle);
    }

    // Ring vertices.
    for (size_t ni = 0; ni < num_nodes; ++ni) {
        const vsg::dvec3 tang = Tangent(pts, ni);
        vsg::dvec3 u, v;
        PerpFrame(tang, u, v);
        for (int s = 0; s < kSides; ++s) {
            const double cs = cos_table[s];
            const double sn = sin_table[s];
            verts[ni * kSides + s] = chrono::ChVector3d(
                pts[ni][0] + kTubeRadius * (cs * u.x + sn * v.x),
                pts[ni][1] + kTubeRadius * (cs * u.y + sn * v.y),
                pts[ni][2] + kTubeRadius * (cs * u.z + sn * v.z));
        }
    }

    // Cap centre vertices.
    verts[ring_verts]     = {pts.front()[0], pts.front()[1], pts.front()[2]};
    verts[ring_verts + 1] = {pts.back()[0],  pts.back()[1],  pts.back()[2]};

    // Tube body quads (two triangles each).
    const size_t body_faces = 2 * (num_nodes - 1) * kSides;
    const size_t cap_faces  = 2 * kSides;
    mesh_faces.reserve(body_faces + cap_faces);

    for (size_t ni = 0; ni < num_nodes - 1; ++ni) {
        for (int s = 0; s < kSides; ++s) {
            const int sn = (s + 1) % kSides;
            const int a  = static_cast<int>(ni       * kSides + s);
            const int b  = static_cast<int>(ni       * kSides + sn);
            const int c  = static_cast<int>((ni + 1) * kSides + sn);
            const int d  = static_cast<int>((ni + 1) * kSides + s);
            mesh_faces.emplace_back(a, d, c);
            mesh_faces.emplace_back(a, c, b);
        }
    }

    // Start cap (fan around ring_verts).
    for (int s = 0; s < kSides; ++s) {
        const int sn = (s + 1) % kSides;
        mesh_faces.emplace_back(static_cast<int>(ring_verts), sn, s);
    }

    // End cap (fan around ring_verts + 1).
    const int end_base = static_cast<int>((num_nodes - 1) * kSides);
    for (int s = 0; s < kSides; ++s) {
        const int sn = (s + 1) % kSides;
        mesh_faces.emplace_back(static_cast<int>(ring_verts + 1),
                                end_base + s, end_base + sn);
    }

    // CreateTrimeshColAvgShape produces indexed geometry with a viewer-
    // connected compile traversal, so subsequent dirty() calls on the vertex
    // buffer trigger GPU re-upload.
    static const chrono::ChColor kCableColor(0.55f, 0.45f, 0.35f);
    auto transform = vsg::MatrixTransform::create();
    lm.node = shape_builder->CreateTrimeshColAvgShape(
        mesh, transform, kCableColor, /*wireframe=*/false);
    if (!lm.node) return;

    // Extract the VSG vertex / normal / colour arrays for per-frame writes.
    lm.vertices = vsg::visit<FindVec3Buffer<0>>(lm.node).Result();
    lm.normals  = vsg::visit<FindVec3Buffer<1>>(lm.node).Result();
    lm.colors   = vsg::visit<FindVec4Buffer<3>>(lm.node).Result();
    if (lm.vertices)
        lm.vertices->properties.dataVariance = vsg::DataVariance::DYNAMIC_DATA;
    if (lm.normals)
        lm.normals->properties.dataVariance = vsg::DataVariance::DYNAMIC_DATA;
    if (lm.colors)
        lm.colors->properties.dataVariance = vsg::DataVariance::DYNAMIC_DATA;

    scene_->addChild(lm.node);

    // Endpoint spheres (colour-coded by MoorDyn point type).
    const auto& p0 = pts.front();
    const auto& pN = pts.back();
    lm.start_marker = CreateSphereMarker(
        vis, p0[0], p0[1], p0[2], kEndpointRadius,
        static_cast<int>(line_data.start_point_type));
    lm.end_marker = CreateSphereMarker(
        vis, pN[0], pN[1], pN[2], kEndpointRadius,
        static_cast<int>(line_data.end_point_type));

    // Small intermediate-node markers to visualise segment boundaries.
    if (num_nodes > 2) {
        lm.node_markers.reserve(num_nodes - 2);
        for (size_t ni = 1; ni + 1 < num_nodes; ++ni) {
            lm.node_markers.push_back(CreateSphereMarker(
                vis, pts[ni][0], pts[ni][1], pts[ni][2],
                kNodeMarkerRadius, kIntermediateNode));
        }
    }
}

// ---------------------------------------------------------------------------

MooringLinesViz::PointMarker MooringLinesViz::CreateSphereMarker(
    chrono::vsg3d::ChVisualSystemVSG* vis,
    double x, double y, double z,
    double radius, int point_type) {
    PointMarker pm;
    auto shape_builder = vis->GetVSGShapeBuilder();
    if (!shape_builder) return pm;

    auto mat = std::make_shared<chrono::ChVisualMaterial>();
    mat->SetDiffuseColor(PointTypeColor(point_type));
    if (point_type == kIntermediateNode) {
        mat->SetMetallic(0.4f);
        mat->SetRoughness(0.4f);
    } else {
        mat->SetMetallic(0.5f);
        mat->SetRoughness(0.3f);
    }

    pm.transform = vsg::MatrixTransform::create();
    pm.transform->matrix =
        vsg::translate(x, y, z) * vsg::scale(radius, radius, radius);

    using ST = chrono::vsg3d::ShapeBuilder::ShapeType;
    pm.node = shape_builder->CreatePbrShape(
        ST::SPHERE, mat, pm.transform, /*wireframe=*/false);
    if (pm.node) scene_->addChild(pm.node);
    return pm;
}

void MooringLinesViz::UpdateMarkerPosition(PointMarker& marker,
                                           double x, double y, double z,
                                           double radius) {
    if (!marker.transform) return;
    marker.transform->matrix =
        vsg::translate(x, y, z) * vsg::scale(radius, radius, radius);
}

// ---------------------------------------------------------------------------

void MooringLinesViz::UpdateTubeMesh(const MooringLineVizData& line_data,
                                     LineMesh& lm,
                                     bool color_enabled,
                                     float range_min, float range_max) {
    const auto& pts        = line_data.node_positions;
    const size_t num_nodes = pts.size();
    if (num_nodes != lm.num_nodes || !lm.vertices) return;

    const size_t ring_verts = num_nodes * static_cast<size_t>(kSides);
    if (lm.vertices->size() < ring_verts + 2) return;

    const bool has_scalars =
        color_enabled && lm.colors &&
        line_data.node_tensions.size() == num_nodes;
    const float inv_range =
        (range_max - range_min > 1e-12f) ? 1.0f / (range_max - range_min)
                                         : 1.0f;

    std::vector<double> cos_table(kSides);
    std::vector<double> sin_table(kSides);
    for (int s = 0; s < kSides; ++s) {
        const double angle = 2.0 * M_PI * s / kSides;
        cos_table[s] = std::cos(angle);
        sin_table[s] = std::sin(angle);
    }

    for (size_t ni = 0; ni < num_nodes; ++ni) {
        const vsg::dvec3 tang = Tangent(pts, ni);
        vsg::dvec3 u, v;
        PerpFrame(tang, u, v);

        const auto centre_x = static_cast<float>(pts[ni][0]);
        const auto centre_y = static_cast<float>(pts[ni][1]);
        const auto centre_z = static_cast<float>(pts[ni][2]);

        vsg::vec4 col = kDefaultCableColorVec4;
        if (has_scalars) {
            const float t =
                (static_cast<float>(line_data.node_tensions[ni]) - range_min)
                * inv_range;
            col = TurboColormap(t);
        }

        for (int s = 0; s < kSides; ++s) {
            const double cs = cos_table[s];
            const double sn = sin_table[s];
            const size_t idx = ni * kSides + s;

            (*lm.vertices)[idx] = vsg::vec3(
                static_cast<float>(pts[ni][0] + kTubeRadius * (cs * u.x + sn * v.x)),
                static_cast<float>(pts[ni][1] + kTubeRadius * (cs * u.y + sn * v.y)),
                static_cast<float>(pts[ni][2] + kTubeRadius * (cs * u.z + sn * v.z)));

            float nx = (*lm.vertices)[idx].x - centre_x;
            float ny = (*lm.vertices)[idx].y - centre_y;
            float nz = (*lm.vertices)[idx].z - centre_z;
            const float nlen = std::sqrt(nx * nx + ny * ny + nz * nz);
            if (nlen > 1e-6f) { nx /= nlen; ny /= nlen; nz /= nlen; }
            (*lm.normals)[idx] = vsg::vec3(nx, ny, nz);

            if (lm.colors && idx < lm.colors->size())
                (*lm.colors)[idx] = col;
        }
    }

    // Cap centres.
    (*lm.vertices)[ring_verts] = vsg::vec3(
        static_cast<float>(pts.front()[0]),
        static_cast<float>(pts.front()[1]),
        static_cast<float>(pts.front()[2]));
    (*lm.vertices)[ring_verts + 1] = vsg::vec3(
        static_cast<float>(pts.back()[0]),
        static_cast<float>(pts.back()[1]),
        static_cast<float>(pts.back()[2]));

    const vsg::dvec3 t0 = Tangent(pts, 0);
    const vsg::dvec3 tN = Tangent(pts, num_nodes - 1);
    (*lm.normals)[ring_verts] = vsg::vec3(
        static_cast<float>(-t0.x), static_cast<float>(-t0.y),
        static_cast<float>(-t0.z));
    (*lm.normals)[ring_verts + 1] = vsg::vec3(
        static_cast<float>(tN.x), static_cast<float>(tN.y),
        static_cast<float>(tN.z));

    if (lm.colors) {
        vsg::vec4 cap0_col = kDefaultCableColorVec4;
        vsg::vec4 capN_col = kDefaultCableColorVec4;
        if (has_scalars) {
            cap0_col = TurboColormap(
                (static_cast<float>(line_data.node_tensions.front()) - range_min)
                * inv_range);
            capN_col = TurboColormap(
                (static_cast<float>(line_data.node_tensions.back()) - range_min)
                * inv_range);
        }
        if (ring_verts < lm.colors->size())
            (*lm.colors)[ring_verts] = cap0_col;
        if (ring_verts + 1 < lm.colors->size())
            (*lm.colors)[ring_verts + 1] = capN_col;
        lm.colors->dirty();
    }

    lm.vertices->dirty();
    lm.normals->dirty();

    // Reposition sphere markers.
    UpdateMarkerPosition(lm.start_marker,
                         pts.front()[0], pts.front()[1], pts.front()[2],
                         kEndpointRadius);
    UpdateMarkerPosition(lm.end_marker,
                         pts.back()[0], pts.back()[1], pts.back()[2],
                         kEndpointRadius);
    for (size_t mi = 0; mi < lm.node_markers.size(); ++mi) {
        const size_t ni = mi + 1;
        if (ni < num_nodes) {
            UpdateMarkerPosition(lm.node_markers[mi],
                                 pts[ni][0], pts[ni][1], pts[ni][2],
                                 kNodeMarkerRadius);
        }
    }
}

}  // namespace gui
}  // namespace hydroc
