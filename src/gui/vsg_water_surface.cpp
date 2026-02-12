// HydroChrono VSG Water Surface Visualization
#include "vsg_water_surface.h"

#include "vsg_config.h"
#include "vsg_gui_component.h"
#include "vsg_materials.h"
#include "vsg_radiation_surface.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>

#include <chrono/core/ChTypes.h>
#include <chrono/physics/ChBodyEasy.h>
#include <hydroc/waves/irregular_wave.h>

#include <Eigen/Dense>

namespace hydroc {
namespace gui {

using namespace vsg_config;

// Visitor to find vertex (N=0) or normal (N=1) arrays in VSG scene graph.
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

// Visitor to find vec4 color array at index N in VSG scene graph.
template <int N>
class FindVec4BufferData : public vsg::Visitor {
  public:
    FindVec4BufferData() : buffer_(nullptr) {}

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

    void apply(vsg::vec4Array& colors) override {
        if (!buffer_) {
            buffer_ = &colors;
        }
    }

    vsg::ref_ptr<vsg::vec4Array> GetBufferData() { return buffer_; }

  private:
    vsg::ref_ptr<vsg::vec4Array> buffer_;
};

// Visitor to enable alpha blending on a VSG graphics pipeline.
// This modifies the ColorBlendState to support transparent rendering.
class EnableAlphaBlending : public vsg::Visitor {
  public:
    EnableAlphaBlending() = default;

    void apply(vsg::Object& object) override { object.traverse(*this); }

    void apply(vsg::GraphicsPipeline& pipeline) override {
        // Find and modify the ColorBlendState in the pipeline's state objects.
        for (auto& state : pipeline.pipelineStates) {
            if (auto* cbs = dynamic_cast<vsg::ColorBlendState*>(state.get())) {
                // Enable alpha blending for the first attachment.
                if (!cbs->attachments.empty()) {
                    auto& att = cbs->attachments[0];
                    att.blendEnable = VK_TRUE;
                    att.srcColorBlendFactor = VK_BLEND_FACTOR_SRC_ALPHA;
                    att.dstColorBlendFactor = VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA;
                    att.colorBlendOp = VK_BLEND_OP_ADD;
                    att.srcAlphaBlendFactor = VK_BLEND_FACTOR_SRC_ALPHA;
                    att.dstAlphaBlendFactor = VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA;
                    att.alphaBlendOp = VK_BLEND_OP_ADD;
                    modified_ = true;
                }
            }
        }
        // Continue traversal in case there are nested pipelines.
        pipeline.traverse(*this);
    }

    bool WasModified() const { return modified_; }

  private:
    bool modified_ = false;
};

AnimatedWaterSurface::AnimatedWaterSurface() = default;
AnimatedWaterSurface::~AnimatedWaterSurface() = default;

void AnimatedWaterSurface::Reset() {
    // Remove water surface from scene if still attached.
    if (scene_ && vsg_node_) {
        auto& children = scene_->children;
        auto it = std::find(children.begin(), children.end(), vsg_node_);
        if (it != children.end()) {
            children.erase(it);
        }
    }

    // Remove wireframe from scene if attached.
    if (scene_ && wireframe_node_) {
        auto& children = scene_->children;
        auto it = std::find(children.begin(), children.end(), wireframe_node_);
        if (it != children.end()) {
            children.erase(it);
        }
    }

    mesh_.reset();
    vsg_vertices_.reset();
    vsg_normals_.reset();
    vsg_colors_.reset();
    vsg_node_.reset();

    wireframe_vertices_.reset();
    wireframe_node_.reset();
    wireframe_visible_ = false;
    wireframe_initialized_ = false;

    scene_.reset();
    bound_vis_ = nullptr;
    initialized_ = false;
    visible_ = true;
    num_triangles_ = 0;
    current_resolution_ = 0;
    last_update_time_ = -1.0;

    // Reset adaptive height-shading range.
    adaptive_eta_min_ = 0.0f;
    adaptive_eta_max_ = 0.0f;
    adaptive_range_initialized_ = false;
}

void AnimatedWaterSurface::Initialize(chrono::vsg3d::ChVisualSystemVSG* vis, int resolution,
                                       const ViewerSettings* settings) {
    // Use default resolution from config if not specified.
    int res = (resolution > 0) ? resolution : kWaterGridResolution;
    InitializeInternal(vis, res, settings);
}

void AnimatedWaterSurface::Reinitialize(int resolution, const ViewerSettings* settings) {
    if (!bound_vis_) {
        return;
    }

    // Check if extent changed.
    bool extent_changed = false;
    if (settings) {
        if (settings->grid_width > 0 && settings->grid_width != grid_width_) extent_changed = true;
        if (settings->grid_length > 0 && settings->grid_length != grid_length_) extent_changed = true;
        if (settings->grid_center_x != grid_center_x_) extent_changed = true;
        if (settings->grid_center_y != grid_center_y_) extent_changed = true;
    }

    // Skip if resolution and extent are the same.
    if (resolution == current_resolution_ && initialized_ && !extent_changed) {
        return;
    }

    if constexpr (kDebugWaveSurface) {
        std::cout << "[WaveSurfaceDebug] Reinitialize with resolution " << resolution << std::endl;
    }

    // Save visual system pointer before reset.
    auto* vis = bound_vis_;

    // Reset and reinitialize with new resolution.
    Reset();
    InitializeInternal(vis, resolution, settings);
}

void AnimatedWaterSurface::InitializeInternal(chrono::vsg3d::ChVisualSystemVSG* vis, int resolution,
                                               const ViewerSettings* settings) {
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

    // Store custom grid extent from settings.
    if (settings) {
        grid_width_ = settings->grid_width;
        grid_length_ = settings->grid_length;
        grid_center_x_ = settings->grid_center_x;
        grid_center_y_ = settings->grid_center_y;
    }

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

    // Use custom grid extent if specified, otherwise use defaults.
    const double grid_size_x = (grid_width_ > 0) ? grid_width_ : kWaterGridSize;
    const double grid_size_y = (grid_length_ > 0) ? grid_length_ : kWaterGridSize;
    const double half_size_x = grid_size_x / 2.0;
    const double half_size_y = grid_size_y / 2.0;
    const double spacing_x = grid_size_x / (n - 1);
    const double spacing_y = grid_size_y / (n - 1);

    // Create vertices (n x n grid centered at custom center).
    std::vector<chrono::ChVector3d>& verts = mesh_->m_vertices;
    verts.resize(static_cast<size_t>(n * n));
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            double x = grid_center_x_ - half_size_x + i * spacing_x;
            double y = grid_center_y_ - half_size_y + j * spacing_y;
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

    // Create transparent white material for the water surface.
    // White diffuse allows vertex colors to show through when multiplied.
    // Opacity enables transparency so we can see through the water.
    auto water_material = chrono_types::make_shared<chrono::ChVisualMaterial>();
    water_material->SetDiffuseColor(chrono::ChColor(1.0f, 1.0f, 1.0f));
    water_material->SetOpacity(kWaterOpacity);
    water_material->SetRoughness(kWaterRoughness);
    water_material->SetMetallic(kWaterMetallic);

    auto transform = vsg::MatrixTransform::create();
    std::vector<chrono::ChVisualMaterialSharedPtr> materials = {water_material};
    vsg_node_ = shape_builder->CreateTrimeshPbrMatShape(mesh_, transform, materials, false);

    if (!vsg_node_) {
        if constexpr (kDebugWaveSurface) {
            std::cout << "[WaveSurfaceDebug] Failed to create VSG mesh\n";
        }
        return;
    }

    // Use VSG visitor pattern to find vertex and normal arrays.
    // For CreateTrimeshPbrMatShape: [0]=vertices, [1]=normals, [2]=texcoords
    vsg_vertices_ = vsg::visit<FindVec3BufferData<0>>(vsg_node_).GetBufferData();
    vsg_normals_ = vsg::visit<FindVec3BufferData<1>>(vsg_node_).GetBufferData();
    // PbrMatShape may not have vertex colors - try to find them anyway.
    vsg_colors_ = vsg::visit<FindVec4BufferData<3>>(vsg_node_).GetBufferData();

    if (!vsg_vertices_ || !vsg_normals_) {
        if constexpr (kDebugWaveSurface) {
            std::cout << "[WaveSurfaceDebug] Failed to find vertex/normal arrays via visitor\n";
        }
        return;
    }

    // Mark arrays as dynamic for GPU re-upload support.
    vsg_vertices_->properties.dataVariance = vsg::DataVariance::DYNAMIC_DATA;
    vsg_normals_->properties.dataVariance = vsg::DataVariance::DYNAMIC_DATA;
    if (vsg_colors_) {
        vsg_colors_->properties.dataVariance = vsg::DataVariance::DYNAMIC_DATA;
    }

    // Add to scene.
    scene_->addChild(vsg_node_);

    initialized_ = true;
    visible_ = true;

    // Log water surface creation.
    const double extent_x = (grid_width_ > 0) ? grid_width_ : kWaterGridSize;
    const double extent_y = (grid_length_ > 0) ? grid_length_ : kWaterGridSize;

    std::cout << "[WaterSurface] Created: " << resolution << "x" << resolution
              << " grid, " << extent_x << " x " << extent_y << " m";
    if (grid_center_x_ != 0.0 || grid_center_y_ != 0.0) {
        std::cout << ", center: (" << grid_center_x_ << ", " << grid_center_y_ << ")";
    }
    std::cout << std::endl;
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

void AnimatedWaterSurface::SetWireframeVisible(bool visible) {
    if (!initialized_ || !scene_) {
        wireframe_visible_ = visible;
        return;
    }

    // Initialize wireframe on first show if needed.
    if (visible && !wireframe_initialized_) {
        InitializeWireframe();
    }

    if (!wireframe_node_ || !wireframe_initialized_) {
        wireframe_visible_ = visible;
        return;
    }

    if (visible == wireframe_visible_) {
        return;  // No change.
    }

    wireframe_visible_ = visible;

    auto& children = scene_->children;
    if (visible) {
        // Add to scene if not present.
        auto it = std::find(children.begin(), children.end(), wireframe_node_);
        if (it == children.end()) {
            children.push_back(wireframe_node_);
        }
    } else {
        // Remove from scene.
        auto it = std::find(children.begin(), children.end(), wireframe_node_);
        if (it != children.end()) {
            children.erase(it);
        }
    }
}

void AnimatedWaterSurface::InitializeWireframe() {
    if (wireframe_initialized_ || !bound_vis_ || !scene_) {
        return;
    }

    auto shape_builder = bound_vis_->GetVSGShapeBuilder();
    if (!shape_builder) {
        return;
    }

    const int n = current_resolution_;
    if (n < 2) {
        return;
    }

    // Create thin quad segments as "lines" that follow the wave surface.
    auto wireframe_mesh = chrono_types::make_shared<chrono::ChTriangleMeshConnected>();

    // Use custom grid extent if specified, otherwise use defaults.
    const double grid_size_x = (grid_width_ > 0) ? grid_width_ : kWaterGridSize;
    const double grid_size_y = (grid_length_ > 0) ? grid_length_ : kWaterGridSize;
    const double half_size_x = grid_size_x / 2.0;
    const double half_size_y = grid_size_y / 2.0;
    const double spacing_x = grid_size_x / (n - 1);
    const double spacing_y = grid_size_y / (n - 1);
    const double line_half_width = 0.05;

    std::vector<chrono::ChVector3d>& verts = wireframe_mesh->m_vertices;
    std::vector<chrono::ChVector3i>& faces = wireframe_mesh->m_face_v_indices;

    int vertex_idx = 0;

    // Horizontal segments (along X, one per grid edge).
    for (int j = 0; j < n; ++j) {
        double y = grid_center_y_ - half_size_y + j * spacing_y;
        for (int i = 0; i < n - 1; ++i) {
            double x0 = grid_center_x_ - half_size_x + i * spacing_x;
            double x1 = grid_center_x_ - half_size_x + (i + 1) * spacing_x;

            verts.push_back(chrono::ChVector3d(x0, y - line_half_width, 0.0));
            verts.push_back(chrono::ChVector3d(x1, y - line_half_width, 0.0));
            verts.push_back(chrono::ChVector3d(x1, y + line_half_width, 0.0));
            verts.push_back(chrono::ChVector3d(x0, y + line_half_width, 0.0));

            faces.push_back(chrono::ChVector3i(vertex_idx, vertex_idx + 1, vertex_idx + 2));
            faces.push_back(chrono::ChVector3i(vertex_idx, vertex_idx + 2, vertex_idx + 3));
            vertex_idx += 4;
        }
    }

    // Vertical segments (along Y, one per grid edge).
    for (int i = 0; i < n; ++i) {
        double x = grid_center_x_ - half_size_x + i * spacing_x;
        for (int j = 0; j < n - 1; ++j) {
            double y0 = grid_center_y_ - half_size_y + j * spacing_y;
            double y1 = grid_center_y_ - half_size_y + (j + 1) * spacing_y;

            verts.push_back(chrono::ChVector3d(x - line_half_width, y0, 0.0));
            verts.push_back(chrono::ChVector3d(x + line_half_width, y0, 0.0));
            verts.push_back(chrono::ChVector3d(x + line_half_width, y1, 0.0));
            verts.push_back(chrono::ChVector3d(x - line_half_width, y1, 0.0));

            faces.push_back(chrono::ChVector3i(vertex_idx, vertex_idx + 1, vertex_idx + 2));
            faces.push_back(chrono::ChVector3i(vertex_idx, vertex_idx + 2, vertex_idx + 3));
            vertex_idx += 4;
        }
    }

    // Faint blue-gray material.
    auto wire_material = chrono_types::make_shared<chrono::ChVisualMaterial>();
    wire_material->SetDiffuseColor(chrono::ChColor(0.1f, 0.2f, 0.3f));
    wire_material->SetOpacity(0.35f);
    wire_material->SetRoughness(0.9f);
    wire_material->SetMetallic(0.0f);

    auto transform = vsg::MatrixTransform::create();
    std::vector<chrono::ChVisualMaterialSharedPtr> materials = {wire_material};
    wireframe_node_ = shape_builder->CreateTrimeshPbrMatShape(wireframe_mesh, transform, materials, false);

    if (!wireframe_node_) {
        return;
    }

    wireframe_vertices_ = vsg::visit<FindVec3BufferData<0>>(wireframe_node_).GetBufferData();
    if (wireframe_vertices_) {
        wireframe_vertices_->properties.dataVariance = vsg::DataVariance::DYNAMIC_DATA;
    }

    wireframe_initialized_ = true;
    std::cout << "[WaterSurface] Wireframe overlay enabled" << std::endl;
}

void AnimatedWaterSurface::UpdateWireframe() {
    if (!wireframe_initialized_ || !wireframe_vertices_ || !vsg_vertices_) {
        return;
    }

    const int n = current_resolution_;
    if (n < 2) {
        return;
    }

    // Get Z at grid position (i, j) from the water surface soup mesh.
    auto get_z_at_grid = [&](int i, int j) -> float {
        i = std::max(0, std::min(i, n - 1));
        j = std::max(0, std::min(j, n - 1));

        // Try as v00 of cell (i, j).
        if (i < n - 1 && j < n - 1) {
            size_t cell_idx = static_cast<size_t>(j * (n - 1) + i);
            size_t idx = cell_idx * 6;
            if (idx < vsg_vertices_->size()) {
                return (*vsg_vertices_)[idx].z;
            }
        }
        // Try as v10 of cell (i-1, j).
        if (i > 0 && j < n - 1) {
            size_t cell_idx = static_cast<size_t>(j * (n - 1) + (i - 1));
            size_t idx = cell_idx * 6 + 1;
            if (idx < vsg_vertices_->size()) {
                return (*vsg_vertices_)[idx].z;
            }
        }
        // Try as v01 of cell (i, j-1).
        if (i < n - 1 && j > 0) {
            size_t cell_idx = static_cast<size_t>((j - 1) * (n - 1) + i);
            size_t idx = cell_idx * 6 + 5;
            if (idx < vsg_vertices_->size()) {
                return (*vsg_vertices_)[idx].z;
            }
        }
        // Try as v11 of cell (i-1, j-1).
        if (i > 0 && j > 0) {
            size_t cell_idx = static_cast<size_t>((j - 1) * (n - 1) + (i - 1));
            size_t idx = cell_idx * 6 + 2;
            if (idx < vsg_vertices_->size()) {
                return (*vsg_vertices_)[idx].z;
            }
        }
        return 0.0f;
    };

    // Use custom grid extent if specified, otherwise use defaults.
    const double grid_size_x = (grid_width_ > 0) ? grid_width_ : kWaterGridSize;
    const double grid_size_y = (grid_length_ > 0) ? grid_length_ : kWaterGridSize;
    const double half_size_x = grid_size_x / 2.0;
    const double half_size_y = grid_size_y / 2.0;
    const double spacing_x = grid_size_x / (n - 1);
    const double spacing_y = grid_size_y / (n - 1);
    const double line_half_width = 0.05;
    const float z_offset = 0.02f;

    size_t soup_idx = 0;

    // Update horizontal segments.
    for (int j = 0; j < n; ++j) {
        double y = grid_center_y_ - half_size_y + j * spacing_y;
        for (int i = 0; i < n - 1; ++i) {
            double x0 = grid_center_x_ - half_size_x + i * spacing_x;
            double x1 = grid_center_x_ - half_size_x + (i + 1) * spacing_x;
            float z0 = get_z_at_grid(i, j) + z_offset;
            float z1 = get_z_at_grid(i + 1, j) + z_offset;

            if (soup_idx + 5 < wireframe_vertices_->size()) {
                (*wireframe_vertices_)[soup_idx + 0] = vsg::vec3(static_cast<float>(x0), static_cast<float>(y - line_half_width), z0);
                (*wireframe_vertices_)[soup_idx + 1] = vsg::vec3(static_cast<float>(x1), static_cast<float>(y - line_half_width), z1);
                (*wireframe_vertices_)[soup_idx + 2] = vsg::vec3(static_cast<float>(x1), static_cast<float>(y + line_half_width), z1);
                (*wireframe_vertices_)[soup_idx + 3] = vsg::vec3(static_cast<float>(x0), static_cast<float>(y - line_half_width), z0);
                (*wireframe_vertices_)[soup_idx + 4] = vsg::vec3(static_cast<float>(x1), static_cast<float>(y + line_half_width), z1);
                (*wireframe_vertices_)[soup_idx + 5] = vsg::vec3(static_cast<float>(x0), static_cast<float>(y + line_half_width), z0);
            }
            soup_idx += 6;
        }
    }

    // Update vertical segments.
    for (int i = 0; i < n; ++i) {
        double x = grid_center_x_ - half_size_x + i * spacing_x;
        for (int j = 0; j < n - 1; ++j) {
            double y0 = grid_center_y_ - half_size_y + j * spacing_y;
            double y1 = grid_center_y_ - half_size_y + (j + 1) * spacing_y;
            float z0 = get_z_at_grid(i, j) + z_offset;
            float z1 = get_z_at_grid(i, j + 1) + z_offset;

            if (soup_idx + 5 < wireframe_vertices_->size()) {
                (*wireframe_vertices_)[soup_idx + 0] = vsg::vec3(static_cast<float>(x - line_half_width), static_cast<float>(y0), z0);
                (*wireframe_vertices_)[soup_idx + 1] = vsg::vec3(static_cast<float>(x + line_half_width), static_cast<float>(y0), z0);
                (*wireframe_vertices_)[soup_idx + 2] = vsg::vec3(static_cast<float>(x + line_half_width), static_cast<float>(y1), z1);
                (*wireframe_vertices_)[soup_idx + 3] = vsg::vec3(static_cast<float>(x - line_half_width), static_cast<float>(y0), z0);
                (*wireframe_vertices_)[soup_idx + 4] = vsg::vec3(static_cast<float>(x + line_half_width), static_cast<float>(y1), z1);
                (*wireframe_vertices_)[soup_idx + 5] = vsg::vec3(static_cast<float>(x - line_half_width), static_cast<float>(y1), z1);
            }
            soup_idx += 6;
        }
    }

    wireframe_vertices_->dirty();
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

    // Check if radiation visualization is enabled.
    // Note: Params are set in UpdateRadiationSourceBody() BEFORE SetSourceState()
    // to ensure amplitude computation uses the correct source radius.
    const bool use_radiation = settings && settings->show_radiation_viz;

    const auto& orig_verts = mesh_->m_vertices;
    const auto& faces = mesh_->m_face_v_indices;

    // Track min/max eta this frame for adaptive range.
    float frame_min_eta = std::numeric_limits<float>::max();
    float frame_max_eta = std::numeric_limits<float>::lowest();

    // Compute adaptive range for height-based shading.
    // Uses previous frame's smoothed range (1-frame lag is imperceptible).
    constexpr float kMinRange = 0.05f;  // Minimum range to avoid extreme sensitivity
    float adaptive_range = std::max(kMinRange, adaptive_eta_max_ - adaptive_eta_min_);
    float adaptive_center = (adaptive_eta_max_ + adaptive_eta_min_) * 0.5f;

    // Debug statistics.
    [[maybe_unused]] double min_eta = std::numeric_limits<double>::max();
    [[maybe_unused]] double max_eta = std::numeric_limits<double>::lowest();

    // -------------------------------------------------------------------------
    // Step 1: Compute wave elevation η(x,y,t) for each unique grid vertex.
    //
    // Performance optimization: The GPU mesh uses a "triangle soup" format where
    // each triangle stores its own copy of vertices. This means interior grid
    // points appear in up to 6 triangles. By computing η once per unique grid
    // point and then copying to the soup, we avoid redundant wave calculations.
    //
    // For irregular sea states with ~1000 frequency components:
    //   - Without optimization: ~24,000 × 1000 = 24M trig operations/frame
    //   - With optimization:     ~4,000 × 1000 = 4M trig operations/frame
    // -------------------------------------------------------------------------
    std::vector<float> grid_eta(orig_verts.size(), 0.0f);

    // For irregular waves, use limited frequency components for visualization.
    // 500 components provides excellent visual quality while being faster than
    // using all components (which can be 1000+). Physics always uses all components.
    constexpr int kVizFrequencyLimit = 500;
    
    // Try to cast to IrregularWaves for the optimized visualization path.
    IrregularWaves* irreg_wave = nullptr;
    if (wave && wave->GetWaveMode() == WaveMode::irregular) {
        irreg_wave = dynamic_cast<IrregularWaves*>(wave.get());
    }

    for (size_t i = 0; i < orig_verts.size(); ++i) {
        const auto& grid_pos = orig_verts[i];
        const Eigen::Vector3d pos(grid_pos.x(), grid_pos.y(), 0.0);

        // Incident wave elevation from the wave model.
        // For irregular waves, use the visualization-optimized method with limited components.
        double eta_incident = 0.0;
        if (irreg_wave) {
            eta_incident = irreg_wave->GetElevationForVisualization(pos, t, kVizFrequencyLimit);
        } else if (wave) {
            eta_incident = wave->GetElevation(pos, t);
        }

        // Radiated wave visualization (approximate, for visual feedback only).
        double eta_radiation = 0.0;
        if (use_radiation && radiation_viz_.IsInitialized()) {
            eta_radiation = radiation_viz_.EvaluateEta(grid_pos.x(), grid_pos.y(), t);
        }

        // Total surface elevation with visual scaling.
        const double eta = (eta_incident + eta_radiation) * visual_scale;
        const float eta_f = static_cast<float>(eta);

        grid_eta[i] = eta_f;

        // Track min/max for adaptive color shading range.
        frame_min_eta = std::min(frame_min_eta, eta_f);
        frame_max_eta = std::max(frame_max_eta, eta_f);

        if constexpr (kDebugWaveSurface) {
            min_eta = std::min(min_eta, eta);
            max_eta = std::max(max_eta, eta);
        }
    }

    // -------------------------------------------------------------------------
    // Step 2: Copy pre-computed elevations to triangle soup vertex buffer.
    // -------------------------------------------------------------------------
    for (size_t tri = 0; tri < faces.size(); ++tri) {
        const auto& face = faces[tri];

        for (int corner = 0; corner < 3; ++corner) {
            int grid_idx = face[corner];
            if (grid_idx < 0 || static_cast<size_t>(grid_idx) >= orig_verts.size()) {
                continue;
            }

            const auto& grid_pos = orig_verts[static_cast<size_t>(grid_idx)];
            float eta_f = grid_eta[static_cast<size_t>(grid_idx)];

            // Write updated position to VSG vertex buffer.
            size_t vsg_idx = tri * 3 + static_cast<size_t>(corner);
            if (vsg_idx < vsg_vertices_->size()) {
                (*vsg_vertices_)[vsg_idx].x = static_cast<float>(grid_pos.x());
                (*vsg_vertices_)[vsg_idx].y = static_cast<float>(grid_pos.y());
                (*vsg_vertices_)[vsg_idx].z = eta_f;
            }

            // Height-based color shading: lighter for crests, darker for troughs.
            if (vsg_colors_ && vsg_idx < vsg_colors_->size()) {
                float norm_eta = (eta_f - adaptive_center) / (adaptive_range * 0.5f);
                norm_eta = std::clamp(norm_eta, -1.0f, 1.0f);
                
                constexpr float kAddOffset = 0.18f;
                float add_term = kAddOffset * norm_eta;
                float brightness = 1.0f + 0.45f * norm_eta;
                
                float r = std::clamp(kWaterR * brightness + add_term * 0.3f, 0.0f, 1.0f);
                float g = std::clamp(kWaterG * brightness + add_term * 0.85f, 0.0f, 1.0f);
                float b = std::clamp(kWaterB * brightness + add_term * 1.0f, 0.0f, 1.0f);
                
                (*vsg_colors_)[vsg_idx] = vsg::vec4(r, g, b, kWaterOpacity);
            }
        }
    }

    // Update adaptive range using exponential smoothing.
    // Fast adaptation when range expands, slow decay when range shrinks.
    // This prevents flickering while still responding to new wave conditions.
    if (frame_min_eta < std::numeric_limits<float>::max()) {
        constexpr float kExpandAlpha = 0.3f;   // Fast expansion (30% new value)
        constexpr float kShrinkAlpha = 0.02f;  // Slow shrinkage (2% new value)
        
        if (!adaptive_range_initialized_) {
            // Initialize on first valid frame.
            adaptive_eta_min_ = frame_min_eta;
            adaptive_eta_max_ = frame_max_eta;
            adaptive_range_initialized_ = true;
        } else {
            // Expand quickly, shrink slowly.
            float min_alpha = (frame_min_eta < adaptive_eta_min_) ? kExpandAlpha : kShrinkAlpha;
            float max_alpha = (frame_max_eta > adaptive_eta_max_) ? kExpandAlpha : kShrinkAlpha;
            
            adaptive_eta_min_ += min_alpha * (frame_min_eta - adaptive_eta_min_);
            adaptive_eta_max_ += max_alpha * (frame_max_eta - adaptive_eta_max_);
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

    // Mark vertex, normal, and color buffers as dirty for GPU re-upload.
    vsg_vertices_->dirty();
    vsg_normals_->dirty();
    if (vsg_colors_) {
        vsg_colors_->dirty();
    }

    // Handle wireframe visibility toggle.
    if (settings) {
        bool want_wireframe = settings->show_water_grid;
        if (want_wireframe != wireframe_visible_) {
            SetWireframeVisible(want_wireframe);
        }
    }

    // Update wireframe if visible.
    if (wireframe_visible_ && wireframe_initialized_) {
        UpdateWireframe();
    }

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

std::shared_ptr<chrono::ChBody> CreateStaticWaterPlane() {
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

