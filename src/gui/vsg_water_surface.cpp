// =============================================================================
// HydroChrono VSG Water Surface Visualization - Implementation
// =============================================================================
#include "vsg_water_surface.h"

#include "vsg_config.h"
#include "vsg_gui_component.h"
#include "vsg_materials.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>

#include <chrono/core/ChTypes.h>
#include <chrono/physics/ChBodyEasy.h>

#include <Eigen/Dense>

namespace hydroc {
namespace gui {

using namespace vsg_config;

// =============================================================================
// Helper: VSG Buffer Visitor
// =============================================================================
// VSG stores geometry data in a tree structure. This "visitor" walks the tree
// to find the actual float arrays containing vertex positions (N=0) or
// normals (N=1). We need these pointers to update the mesh each frame.
//
// This follows the same pattern used internally by Chrono-VSG.
// =============================================================================
template <int N>
class FindVec3BufferData : public vsg::Visitor {
  public:
    FindVec3BufferData() : buffer_(nullptr) {}

    void apply(vsg::Object& object) override { object.traverse(*this); }

    void apply(vsg::VertexDraw& vd) override {
        if (vd.arrays.empty() || static_cast<int>(vd.arrays.size()) <= N) {
            return;
        }
        vd.arrays[N]->data->accept(*this);
    }

    void apply(vsg::VertexIndexDraw& vid) override {
        if (vid.arrays.empty() || static_cast<int>(vid.arrays.size()) <= N) {
            return;
        }
        vid.arrays[N]->data->accept(*this);
    }

    void apply(vsg::vec3Array& vertices) override {
        if (!buffer_) {
            buffer_ = &vertices;
        }
    }

    vsg::ref_ptr<vsg::vec3Array> GetBufferData() { return buffer_; }

  private:
    vsg::ref_ptr<vsg::vec3Array> buffer_;
};

// =============================================================================
// AnimatedWaterSurface Implementation
// =============================================================================

AnimatedWaterSurface::AnimatedWaterSurface() = default;
AnimatedWaterSurface::~AnimatedWaterSurface() = default;

void AnimatedWaterSurface::Reset() {
    // Remove from scene if still attached.
    if (scene_ && vsg_node_) {
        auto& children = scene_->children;
        auto it = std::find(children.begin(), children.end(), vsg_node_);
        if (it != children.end()) {
            children.erase(it);
        }
    }

    mesh_.reset();
    vsg_vertices_.reset();
    vsg_normals_.reset();
    vsg_node_.reset();
    scene_.reset();
    bound_vis_ = nullptr;
    initialized_ = false;
    visible_ = true;
    num_triangles_ = 0;
    current_resolution_ = 0;
    last_update_time_ = -1.0;
}

void AnimatedWaterSurface::Initialize(chrono::vsg3d::ChVisualSystemVSG* vis, int resolution) {
    // Use default resolution from config if not specified.
    int res = (resolution > 0) ? resolution : kWaterGridResolution;
    InitializeInternal(vis, res);
}

void AnimatedWaterSurface::Reinitialize(int resolution) {
    if (!bound_vis_) {
        return;
    }

    // Skip if resolution is the same.
    if (resolution == current_resolution_ && initialized_) {
        return;
    }

    if constexpr (kDebugWaveSurface) {
        std::cout << "[WaveSurfaceDebug] Reinitialize with resolution " << resolution << std::endl;
    }

    // Save visual system pointer before reset.
    auto* vis = bound_vis_;

    // Reset and reinitialize with new resolution.
    Reset();
    InitializeInternal(vis, resolution);
}

void AnimatedWaterSurface::InitializeInternal(chrono::vsg3d::ChVisualSystemVSG* vis, int resolution) {
    if constexpr (kDebugWaveSurface) {
        std::cout << "[WaveSurfaceDebug] Initialize called, vis=" << static_cast<void*>(vis)
                  << " resolution=" << resolution << std::endl;
    }
    if (!vis) {
        return;
    }

    // If visual system changed, reset and reinitialize.
    if (initialized_ && bound_vis_ != vis) {
        if constexpr (kDebugWaveSurface) {
            std::cout << "[WaveSurfaceDebug] vis changed, resetting" << std::endl;
        }
        Reset();
    }
    if (initialized_) {
        return;
    }

    bound_vis_ = vis;
    current_resolution_ = resolution;

    auto shape_builder = vis->GetVSGShapeBuilder();
    scene_ = vis->GetVSGScene();
    if (!shape_builder || !scene_) {
        if constexpr (kDebugWaveSurface) {
            std::cout << "[WaveSurfaceDebug] Failed to get shapeBuilder or scene\n";
        }
        return;
    }

    // Create the Chrono mesh geometry.
    mesh_ = chrono_types::make_shared<chrono::ChTriangleMeshConnected>();

    const int n = resolution;
    const double half_size = kWaterGridSize / 2.0;
    const double spacing = kWaterGridSize / (n - 1);

    // Create vertices (n x n grid centered at origin).
    std::vector<chrono::ChVector3d>& verts = mesh_->m_vertices;
    verts.resize(static_cast<size_t>(n * n));
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            double x = -half_size + i * spacing;
            double y = -half_size + j * spacing;
            verts[static_cast<size_t>(j * n + i)] = chrono::ChVector3d(x, y, 0.0);
        }
    }

    // Create triangle faces (two triangles per grid cell).
    std::vector<chrono::ChVector3i>& faces = mesh_->m_face_v_indices;
    faces.reserve(static_cast<size_t>(2 * (n - 1) * (n - 1)));
    for (int j = 0; j < n - 1; ++j) {
        for (int i = 0; i < n - 1; ++i) {
            int v00 = j * n + i;
            int v10 = j * n + i + 1;
            int v01 = (j + 1) * n + i;
            int v11 = (j + 1) * n + i + 1;
            // Two triangles per quad.
            faces.push_back(chrono::ChVector3i(v00, v10, v11));
            faces.push_back(chrono::ChVector3i(v00, v11, v01));
        }
    }

    num_triangles_ = faces.size();

    // Create water material with translucency and PBR properties.
    auto water_material = chrono_types::make_shared<chrono::ChVisualMaterial>();
    water_material->SetDiffuseColor(chrono::ChColor(kWaterR, kWaterG, kWaterB));
    water_material->SetOpacity(kWaterOpacity);
    water_material->SetRoughness(kWaterRoughness);
    water_material->SetMetallic(kWaterMetallic);
    water_material->SetSpecularColor(chrono::ChColor(kWaterSpecular, kWaterSpecular, kWaterSpecular));

    // Use CreateTrimeshPbrMatShape for proper material with transparency.
    // This creates dynamic-compatible buffers with PBR rendering.
    auto transform = vsg::MatrixTransform::create();
    std::vector<chrono::ChVisualMaterialSharedPtr> materials = {water_material};
    vsg_node_ = shape_builder->CreateTrimeshPbrMatShape(mesh_, transform, materials, false);

    if (!vsg_node_) {
        if constexpr (kDebugWaveSurface) {
            std::cout << "[WaveSurfaceDebug] Failed to create VSG mesh\n";
        }
        return;
    }

    // Use VSG visitor pattern to find vertex and normal arrays (same as Chrono-VSG).
    vsg_vertices_ = vsg::visit<FindVec3BufferData<0>>(vsg_node_).GetBufferData();
    vsg_normals_ = vsg::visit<FindVec3BufferData<1>>(vsg_node_).GetBufferData();

    if (!vsg_vertices_ || !vsg_normals_) {
        if constexpr (kDebugWaveSurface) {
            std::cout << "[WaveSurfaceDebug] Failed to find vertex/normal arrays via visitor\n";
        }
        return;
    }

    // Mark arrays as dynamic for GPU re-upload support.
    vsg_vertices_->properties.dataVariance = vsg::DataVariance::DYNAMIC_DATA;
    vsg_normals_->properties.dataVariance = vsg::DataVariance::DYNAMIC_DATA;

    // Add to scene.
    scene_->addChild(vsg_node_);

    initialized_ = true;
    visible_ = true;

    // Always print one-time status log when water surface is created.
    std::cout << "[WaterSurface] init vis=" << static_cast<void*>(vis)
              << " verts=" << vsg_vertices_->size()
              << " normals=" << vsg_normals_->size()
              << " triangles=" << num_triangles_
              << " resolution=" << resolution
              << " bound=yes" << std::endl;
}

void AnimatedWaterSurface::SetVisible(bool visible) {
    if (!initialized_ || !vsg_node_ || !scene_) {
        visible_ = visible;
        return;
    }

    if (visible == visible_) {
        return;  // No change.
    }

    visible_ = visible;

    auto& children = scene_->children;
    if (visible) {
        // Add back to scene if not present.
        auto it = std::find(children.begin(), children.end(), vsg_node_);
        if (it == children.end()) {
            children.push_back(vsg_node_);
            if constexpr (kDebugWaveSurface) {
                std::cout << "[WaveSurfaceDebug] Water surface shown\n";
            }
        }
    } else {
        // Remove from scene.
        auto it = std::find(children.begin(), children.end(), vsg_node_);
        if (it != children.end()) {
            children.erase(it);
            if constexpr (kDebugWaveSurface) {
                std::cout << "[WaveSurfaceDebug] Water surface hidden\n";
            }
        }
    }
}

bool AnimatedWaterSurface::IsVisible() const {
    return visible_;
}

void AnimatedWaterSurface::Update(const std::shared_ptr<WaveBase>& wave, double t,
                                  const ViewerSettings* settings) {
    if (!initialized_ || !vsg_vertices_ || !vsg_normals_ || !mesh_) {
        if constexpr (kDebugWaveSurface) {
            if (frame_count_++ % kDebugPrintEveryNFrames == 0) {
                std::cout << "[WaveSurfaceDebug] Update() early return: initialized_="
                          << initialized_ << ", vsg_vertices_=" << (vsg_vertices_ ? "valid" : "null")
                          << "\n";
            }
        }
        return;
    }

    // Skip updates if not visible.
    if (!visible_) {
        return;
    }

    // Throttle updates based on update_hz setting.
    if (settings && settings->update_hz > 0) {
        double min_interval = 1.0 / settings->update_hz;
        if (last_update_time_ >= 0.0 && (t - last_update_time_) < min_interval) {
            return;  // Skip this update (throttled).
        }
    }
    last_update_time_ = t;

    // Get visual scale from settings (default 1.0).
    float visual_scale = (settings) ? settings->wave_visual_scale : 1.0f;

    const auto& orig_verts = mesh_->m_vertices;
    const auto& faces = mesh_->m_face_v_indices;

    // Debug statistics (only computed when debug is enabled).
    [[maybe_unused]] double min_eta = std::numeric_limits<double>::max();
    [[maybe_unused]] double max_eta = std::numeric_limits<double>::lowest();

    // Update each triangle's vertices with wave elevation.
    // "Mesh soup" format: each triangle has 3 separate (non-shared) vertices,
    // laid out as [tri0_v0, tri0_v1, tri0_v2, tri1_v0, tri1_v1, tri1_v2, ...]
    for (size_t tri = 0; tri < faces.size(); ++tri) {
        const auto& face = faces[tri];

        for (int corner = 0; corner < 3; ++corner) {
            int grid_idx = face[corner];  // Index into original grid vertices
            if (grid_idx < 0 || static_cast<size_t>(grid_idx) >= orig_verts.size()) {
                continue;
            }

            const auto& grid_pos = orig_verts[static_cast<size_t>(grid_idx)];

            // Get wave elevation eta(x, y, t) at this grid point.
            double eta = 0.0;
            if (wave) {
                Eigen::Vector3d pos(grid_pos.x(), grid_pos.y(), 0.0);
                eta = wave->GetElevation(pos, t);
            }

            // Apply visual scale multiplier.
            eta *= visual_scale;

            // Write updated position to VSG vertex buffer.
            size_t vsg_idx = tri * 3 + static_cast<size_t>(corner);
            if (vsg_idx < vsg_vertices_->size()) {
                (*vsg_vertices_)[vsg_idx].x = static_cast<float>(grid_pos.x());
                (*vsg_vertices_)[vsg_idx].y = static_cast<float>(grid_pos.y());
                (*vsg_vertices_)[vsg_idx].z = static_cast<float>(eta);
            }

            if constexpr (kDebugWaveSurface) {
                min_eta = std::min(min_eta, eta);
                max_eta = std::max(max_eta, eta);
            }
        }
    }

    // Debug: spike first triangle to verify GPU update pipeline is working.
    if constexpr (kDebugSpikeVertex) {
        if (vsg_vertices_->size() >= 3) {
            (*vsg_vertices_)[0].z = 20.0f;
            (*vsg_vertices_)[1].z = 20.0f;
            (*vsg_vertices_)[2].z = 20.0f;
        }
    }

    // Recompute surface normals for correct lighting.
    // Each triangle gets a flat normal computed from cross product of edges.
    for (size_t tri = 0; tri < faces.size(); ++tri) {
        size_t i0 = tri * 3;
        size_t i1 = tri * 3 + 1;
        size_t i2 = tri * 3 + 2;

        if (i2 >= vsg_vertices_->size()) {
            continue;
        }

        const auto& v0 = (*vsg_vertices_)[i0];
        const auto& v1 = (*vsg_vertices_)[i1];
        const auto& v2 = (*vsg_vertices_)[i2];

        // Edge vectors from v0.
        vsg::vec3 edge1(v1.x - v0.x, v1.y - v0.y, v1.z - v0.z);
        vsg::vec3 edge2(v2.x - v0.x, v2.y - v0.y, v2.z - v0.z);

        // Normal = edge1 × edge2 (cross product).
        vsg::vec3 normal(
            edge1.y * edge2.z - edge1.z * edge2.y,
            edge1.z * edge2.x - edge1.x * edge2.z,
            edge1.x * edge2.y - edge1.y * edge2.x);

        // Normalize to unit length.
        float len = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
        if (len > 1e-6f) {
            normal.x /= len;
            normal.y /= len;
            normal.z /= len;
        } else {
            normal = vsg::vec3(0.0f, 0.0f, 1.0f);  // Default to up if degenerate
        }

        // Assign same normal to all 3 vertices of this triangle.
        if (i0 < vsg_normals_->size()) {
            (*vsg_normals_)[i0] = normal;
        }
        if (i1 < vsg_normals_->size()) {
            (*vsg_normals_)[i1] = normal;
        }
        if (i2 < vsg_normals_->size()) {
            (*vsg_normals_)[i2] = normal;
        }
    }

    // Mark vertex and normal buffers as dirty for GPU re-upload.
    vsg_vertices_->dirty();
    vsg_normals_->dirty();

    // Debug output.
    if constexpr (kDebugWaveSurface) {
        if (frame_count_++ % kDebugPrintEveryNFrames == 0) {
            std::cout << "[WaveSurface] t=" << t
                      << " wave=" << (wave ? "ok" : "null")
                      << " scale=" << visual_scale
                      << " min_eta=" << min_eta
                      << " max_eta=" << max_eta;
            if (wave) {
                const double L = kWaterGridSize;
                Eigen::Vector3d p0(0.0, 0.0, 0.0);
                Eigen::Vector3d p1(L / 4.0, 0.0, 0.0);
                Eigen::Vector3d p2(L / 2.0, 0.0, 0.0);
                double eta0 = wave->GetElevation(p0, t);
                double eta1 = wave->GetElevation(p1, t);
                double eta2 = wave->GetElevation(p2, t);
                std::cout << " eta(0,0)=" << eta0
                          << " eta(L/4,0)=" << eta1
                          << " eta(L/2,0)=" << eta2;
            }
            std::cout << "\n";
        }
    }
}

bool AnimatedWaterSurface::IsInitialized() const {
    return initialized_;
}

bool AnimatedWaterSurface::IsInitializedFor(chrono::vsg3d::ChVisualSystemVSG* vis) const {
    return initialized_ && bound_vis_ == vis;
}

std::string AnimatedWaterSurface::GetStatusString() const {
    std::ostringstream ss;
    ss << "Water: " << (initialized_ ? "ANIM" : "OFF");
    ss << " verts=" << (vsg_vertices_ ? std::to_string(vsg_vertices_->size()) : "null");
    return ss.str();
}

// =============================================================================
// Static Water Plane Factory
// =============================================================================

std::shared_ptr<chrono::ChBody> CreateStaticWaterPlane() {
    // Create a thin box as the water surface (visual only: no collision, fixed).
    auto water = chrono_types::make_shared<chrono::ChBodyEasyBox>(
        kWaterPlaneSize,       // x dimension
        kWaterPlaneSize,       // y dimension
        kWaterPlaneThickness,  // z dimension (thin slab)
        1.0,                   // density (irrelevant for fixed body)
        true,                  // enable visualization
        false                  // no collision geometry
    );

    water->SetName("water_surface");
    water->SetPos(chrono::ChVector3d(0.0, 0.0, kWaterPlaneZ));
    water->SetFixed(true);
    water->EnableCollision(false);

    // Apply water material to all visual shapes.
    const auto water_mat = MakeWaterSurfaceMaterial();
    ApplyMaterialToAllVisualShapes(*water, water_mat);

    return water;
}

}  // namespace gui
}  // namespace hydroc

