#include "guihelper_impl.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <vector>

#include <chrono/physics/ChBody.h>
#include <chrono/physics/ChBodyEasy.h>
#include <chrono/assets/ChVisualShapeTriangleMesh.h>
#include <chrono/geometry/ChTriangleMeshConnected.h>

using namespace hydroc::gui;

using namespace chrono::vsg3d;

// -----------------------------------------------------------------------------

namespace {
// -----------------------------------------------------------------------------
// Appearance tuning (Painted metal base — industrial yellow)
// -----------------------------------------------------------------------------
constexpr float kPaintedMetalR = 0.92f;
constexpr float kPaintedMetalG = 0.72f;
constexpr float kPaintedMetalB = 0.12f;

// Moderate, subtle highlight (painted steel: not mirror, not plastic).
constexpr float kSpecular = 0.25f;
constexpr float kSpecularExponent = 64.0f;

// Metallic/roughness parameters (painted metal: mostly dielectric coating).
constexpr float kRoughness = 0.45f;
constexpr float kMetallic = 0.05f;

// Per-body color variation range [0.92, 1.08] for visual distinction.
constexpr float kColorVariationMin = 0.92f;
constexpr float kColorVariationMax = 1.08f;
constexpr int kColorVariationCycle = 8;  // cycle length for deterministic variation

// -----------------------------------------------------------------------------
// Appearance tuning (Water surface)
// -----------------------------------------------------------------------------
// Translucent ocean blue — tropical/deep water aesthetic.
constexpr float kWaterR = 0.01f;
constexpr float kWaterG = 0.20f;
constexpr float kWaterB = 0.35f;

// Translucent surface (allows seeing through to objects below).
constexpr float kWaterOpacity = 0.55f;  // 0=fully transparent, 1=fully opaque

// Specular reflection intensity (white highlight on water surface).
constexpr float kWaterSpecular = 0.6f;

// Glossy surface: low roughness = more reflective (like calm water).
constexpr float kWaterRoughness = 0.05f;
constexpr float kWaterMetallic = 0.0f;  // Water is non-metallic (dielectric)

// Water plane geometry.
constexpr double kWaterPlaneSize = 200.0;       // meters (large enough to fill most views)
constexpr double kWaterPlaneThickness = 0.02;   // meters (thin slab)
constexpr double kWaterPlaneZ = -0.05;          // below z=0 to avoid z-fighting with craft

// -----------------------------------------------------------------------------
// Animated water surface grid (for wave visualization)
// -----------------------------------------------------------------------------
constexpr double kWaterGridSize = 100.0;        // meters (grid extent in x and y)
constexpr int kWaterGridResolution = 64;        // vertices per side (64×64 for smoother surface)
constexpr double kWaterGridSpacing = kWaterGridSize / (kWaterGridResolution - 1);

// -----------------------------------------------------------------------------
// Debug instrumentation (set to true to diagnose flat free-surface issues)
// -----------------------------------------------------------------------------
constexpr bool kDebugWaveSurface = false;       // OFF by default; set true to diagnose
constexpr int kDebugPrintEveryNFrames = 60;     // Print every N frames to reduce spam
constexpr bool kDebugSpikeVertex = false;       // OFF by default; set true to verify VSG pipeline
constexpr bool kForceWaterSurface = false;      // OFF by default; set true to force water even without waves

// -----------------------------------------------------------------------------
// Camera tuning (deterministic side-on view)
// -----------------------------------------------------------------------------
// Eye on negative Y axis, slightly above, looking at origin. Z-up enforced.
constexpr double kCameraDistance = 60.0;        // meters along -Y axis
constexpr double kCameraHeight = 40.0;          // meters above z=0

// -----------------------------------------------------------------------------
// Lighting tuning (3-light studio setup: key + fill + ambient)
// -----------------------------------------------------------------------------
// Key light (main sun) — primary illumination from upper-front.
constexpr float kKeyIntensity = 1.0f;
constexpr double kKeyAzimuth = 0.4 * chrono::CH_PI;      // ~72° from front
constexpr double kKeyElevation = chrono::CH_PI / 5.0;    // ~36° above horizon

// Fill light (soft opposite) — reduces harsh shadows, lifts dark side.
constexpr float kFillIntensity = 0.35f;
constexpr double kFillAzimuth = kKeyAzimuth + chrono::CH_PI;  // opposite side
constexpr double kFillElevation = chrono::CH_PI / 8.0;        // ~22.5° above horizon

// Note: Chrono-VSG includes a built-in ambient at 10% of key intensity.
// Shadows remain OFF (water plane lacks per-object shadow control).

// -----------------------------------------------------------------------------
// Wireframe overlay tuning
// -----------------------------------------------------------------------------
// Enable wireframe rendering for improved shape readability.
// Default OFF: solid shaded rendering displays materials/colors properly.
constexpr bool kEnableWireframe = false;

/// Creates the base painted-metal material (industrial yellow).
chrono::ChVisualMaterial MakePaintedMetalBase() {
    chrono::ChVisualMaterial mat;

    mat.SetDiffuseColor(chrono::ChColor(kPaintedMetalR, kPaintedMetalG, kPaintedMetalB));
    mat.SetSpecularColor(chrono::ChColor(kSpecular, kSpecular, kSpecular));
    mat.SetSpecularExponent(kSpecularExponent);
    mat.SetRoughness(kRoughness);
    mat.SetMetallic(kMetallic);
    mat.SetOpacity(1.0f);

    return mat;
}

/// Creates a painted-metal variant with subtle brightness variation for body distinction.
/// The variation is deterministic based on the body index, cycling through kColorVariationCycle.
chrono::ChVisualMaterial MakePaintedMetalVariant(int index) {
    // Compute variation factor in [kColorVariationMin, kColorVariationMax].
    // Cycle through values to keep variation subtle and deterministic.
    const float t = static_cast<float>(index % kColorVariationCycle) /
                    static_cast<float>(kColorVariationCycle - 1);
    const float factor = kColorVariationMin + t * (kColorVariationMax - kColorVariationMin);

    // Apply factor to base color, clamping to valid range.
    const float r = std::min(1.0f, kPaintedMetalR * factor);
    const float g = std::min(1.0f, kPaintedMetalG * factor);
    const float b = std::min(1.0f, kPaintedMetalB * factor);

    chrono::ChVisualMaterial mat;
    mat.SetDiffuseColor(chrono::ChColor(r, g, b));
    mat.SetSpecularColor(chrono::ChColor(kSpecular, kSpecular, kSpecular));
    mat.SetSpecularExponent(kSpecularExponent);
    mat.SetRoughness(kRoughness);
    mat.SetMetallic(kMetallic);
    mat.SetOpacity(1.0f);

    return mat;
}

/// Creates the water surface material (translucent ocean blue).
chrono::ChVisualMaterial MakeWaterSurfaceMaterial() {
    chrono::ChVisualMaterial mat;

    mat.SetDiffuseColor(chrono::ChColor(kWaterR, kWaterG, kWaterB));
    mat.SetSpecularColor(chrono::ChColor(kWaterSpecular, kWaterSpecular, kWaterSpecular));
    mat.SetRoughness(kWaterRoughness);
    mat.SetMetallic(kWaterMetallic);
    mat.SetOpacity(kWaterOpacity);

    return mat;
}

void ApplyMaterialToAllVisualShapes(chrono::ChBody& body, const chrono::ChVisualMaterial& mat) {
    auto model = body.GetVisualModel();
    if (!model) {
        return;
    }

    auto mat_ptr = chrono_types::make_shared<chrono::ChVisualMaterial>(mat);

    const unsigned int num_shapes = model->GetNumShapes();
    for (unsigned int i = 0; i < num_shapes; ++i) {
        auto shape = model->GetShape(i);
        if (!shape) {
            continue;
        }

        const unsigned int num_materials = shape->GetNumMaterials();
        if (num_materials == 0) {
            shape->AddMaterial(mat_ptr);
            continue;
        }

        for (unsigned int m = 0; m < num_materials; ++m) {
            shape->SetMaterial(m, mat_ptr);
        }
    }
}

/// Adds a fill light to the VSG scene for improved shape readability.
/// Must be called AFTER pVis->Initialize() since the scene is created during init.
void AddFillLight(chrono::vsg3d::ChVisualSystemVSG* vis) {
    auto scene = vis->GetVSGScene();
    if (!scene) {
        return;
    }

    // Compute fill light direction (Z-up coordinate system).
    const double se = std::sin(kFillElevation);
    const double ce = std::cos(kFillElevation);
    const double sa = std::sin(kFillAzimuth);
    const double ca = std::cos(kFillAzimuth);

    auto fillLight = vsg::DirectionalLight::create();
    fillLight->name = "fill light";
    fillLight->color.set(1.0f, 1.0f, 1.0f);
    fillLight->intensity = kFillIntensity;
    fillLight->direction.set(
        static_cast<float>(-ce * ca),
        static_cast<float>(-ce * sa),
        static_cast<float>(-se));

    scene->addChild(fillLight);
}

/// Creates a visual-only water plane (fixed, non-colliding) at z ~ 0.
/// Used when no wave model is set (static still water).
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

}  // namespace

namespace hydroc {
namespace gui {

// =============================================================================
// Helper class to find GPU vertex/normal buffers in VSG scene graph.
// 
// VSG stores geometry data in a tree structure. This "visitor" walks the tree
// to find the actual float arrays containing vertex positions (N=0) or
// normals (N=1). We need these pointers to update the mesh each frame.
//
// This follows the same pattern used internally by Chrono-VSG.
// =============================================================================
template <int N>
class FindVec3BufferData : public vsg::Visitor {
  public:
    FindVec3BufferData() : m_buffer(nullptr) {}
    
    void apply(vsg::Object& object) override { object.traverse(*this); }
    
    void apply(vsg::VertexDraw& vd) override {
        if (vd.arrays.empty() || static_cast<int>(vd.arrays.size()) <= N)
            return;
        vd.arrays[N]->data->accept(*this);
    }
    
    void apply(vsg::VertexIndexDraw& vid) override {
        if (vid.arrays.empty() || static_cast<int>(vid.arrays.size()) <= N)
            return;
        vid.arrays[N]->data->accept(*this);
    }
    
    void apply(vsg::vec3Array& vertices) override {
        if (!m_buffer)
            m_buffer = &vertices;
    }
    
    vsg::ref_ptr<vsg::vec3Array> getBufferData() { return m_buffer; }
    
  private:
    vsg::ref_ptr<vsg::vec3Array> m_buffer;
};

// =============================================================================
// AnimatedWaterSurface
// 
// Manages a deformable water surface mesh for visualizing ocean waves.
// 
// How it works:
// 1. Creates a flat grid of triangles (64×64 vertices by default).
// 2. Each frame, queries the wave model for elevation eta(x,y,t) at each vertex.
// 3. Updates the Z-coordinate of each vertex to match the wave height.
// 4. Recomputes surface normals for correct lighting.
// 5. Tells the GPU to re-upload the modified vertex data.
//
// The mesh uses "soup" format where each triangle has its own 3 vertices
// (not shared with neighbors), which simplifies GPU buffer updates.
// =============================================================================
class AnimatedWaterSurface {
  public:
    AnimatedWaterSurface() = default;

    /// Reset state for reuse with a new system.
    void Reset() {
        mesh_.reset();
        vsg_vertices_.reset();
        vsg_normals_.reset();
        vsg_node_.reset();
        bound_vis_ = nullptr;
        initialized_ = false;
        num_triangles_ = 0;
    }

    /// Initialize the mesh geometry and add to the visual system scene.
    /// Must be called after pVis->Initialize() so the scene exists.
    void Initialize(chrono::vsg3d::ChVisualSystemVSG* vis) {
        if constexpr (kDebugWaveSurface) {
            std::cout << "[WaveSurfaceDebug] Initialize called, vis=" << static_cast<void*>(vis) << std::endl;
        }
        if (!vis) return;

        // If visual system changed, reset and reinitialize.
        if (initialized_ && bound_vis_ != vis) {
            if constexpr (kDebugWaveSurface) {
                std::cout << "[WaveSurfaceDebug] vis changed, resetting" << std::endl;
            }
            Reset();
        }
        if (initialized_) return;

        bound_vis_ = vis;

        auto shapeBuilder = vis->GetVSGShapeBuilder();
        auto scene = vis->GetVSGScene();
        if (!shapeBuilder || !scene) {
            if constexpr (kDebugWaveSurface) {
                std::cout << "[WaveSurfaceDebug] Failed to get shapeBuilder or scene\n";
            }
            return;
        }

        // Create the Chrono mesh geometry.
        mesh_ = chrono_types::make_shared<chrono::ChTriangleMeshConnected>();
        
        const int n = kWaterGridResolution;
        const double half_size = kWaterGridSize / 2.0;

        // Create vertices (n x n grid centered at origin).
        std::vector<chrono::ChVector3d>& verts = mesh_->m_vertices;
        verts.resize(static_cast<size_t>(n * n));
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                double x = -half_size + i * kWaterGridSpacing;
                double y = -half_size + j * kWaterGridSpacing;
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
        vsg_node_ = shapeBuilder->CreateTrimeshPbrMatShape(mesh_, transform, materials, false);

        if (!vsg_node_) {
            if constexpr (kDebugWaveSurface) {
                std::cout << "[WaveSurfaceDebug] Failed to create VSG mesh\n";
            }
            return;
        }

        // Use VSG visitor pattern to find vertex and normal arrays (same as Chrono-VSG).
        vsg_vertices_ = vsg::visit<FindVec3BufferData<0>>(vsg_node_).getBufferData();
        vsg_normals_ = vsg::visit<FindVec3BufferData<1>>(vsg_node_).getBufferData();

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
        scene->addChild(vsg_node_);

        initialized_ = true;

        // Always print one-time status log when water surface is created.
        std::cout << "[WaterSurface] init vis=" << static_cast<void*>(vis)
                  << " verts=" << vsg_vertices_->size()
                  << " normals=" << vsg_normals_->size()
                  << " triangles=" << num_triangles_
                  << " bound=yes" << std::endl;
    }

    /// Update vertex Z coordinates based on wave elevation at time t.
    /// If wave is null, keeps surface flat (useful for kForceWaterSurface mode).
    /// Note: mesh soup format has 3 vertices per triangle (not shared).
    void Update(const std::shared_ptr<WaveBase>& wave, double t) {
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
                if (grid_idx < 0 || static_cast<size_t>(grid_idx) >= orig_verts.size())
                    continue;

                const auto& grid_pos = orig_verts[static_cast<size_t>(grid_idx)];
                
                // Get wave elevation eta(x, y, t) at this grid point.
                double eta = 0.0;
                if (wave) {
                    Eigen::Vector3d pos(grid_pos.x(), grid_pos.y(), 0.0);
                    eta = wave->GetElevation(pos, t);
                }

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
            
            if (i2 >= vsg_vertices_->size()) continue;

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
                edge1.x * edge2.y - edge1.y * edge2.x
            );
            
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
            if (i0 < vsg_normals_->size()) (*vsg_normals_)[i0] = normal;
            if (i1 < vsg_normals_->size()) (*vsg_normals_)[i1] = normal;
            if (i2 < vsg_normals_->size()) (*vsg_normals_)[i2] = normal;
        }

        // Mark vertex and normal buffers as dirty for GPU re-upload.
        vsg_vertices_->dirty();
        vsg_normals_->dirty();

        // Debug output
        if constexpr (kDebugWaveSurface) {
            if (frame_count_++ % kDebugPrintEveryNFrames == 0) {
                std::cout << "[WaveSurface] t=" << t 
                          << " wave=" << (wave ? "ok" : "null")
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

    bool IsInitialized() const { return initialized_; }

    /// Check if initialized for the given visual system.
    bool IsInitializedFor(chrono::vsg3d::ChVisualSystemVSG* vis) const {
        return initialized_ && bound_vis_ == vis;
    }

    /// Get a human-readable status string (for logging).
    std::string GetStatusString() const {
        std::ostringstream ss;
        ss << "Water: " << (initialized_ ? "ANIM" : "OFF");
        ss << " verts=" << (vsg_vertices_ ? std::to_string(vsg_vertices_->size()) : "null");
        return ss.str();
    }

  private:
    std::shared_ptr<chrono::ChTriangleMeshConnected> mesh_;
    vsg::ref_ptr<vsg::vec3Array> vsg_vertices_;   ///< VSG vertex positions (dynamic, mesh soup)
    vsg::ref_ptr<vsg::vec3Array> vsg_normals_;    ///< VSG vertex normals (dynamic, mesh soup)
    vsg::ref_ptr<vsg::Node> vsg_node_;            ///< Reference to VSG scene node
    chrono::vsg3d::ChVisualSystemVSG* bound_vis_ = nullptr;
    size_t num_triangles_ = 0;                    ///< Number of triangles in mesh
    bool initialized_ = false;
    mutable int frame_count_ = 0;  // Frame counter for debug output throttling
};

}  // namespace gui
}  // namespace hydroc

// -----------------------------------------------------------------------------

class MyComponentVSG : public chrono::vsg3d::ChGuiComponentVSG {
  public:
    MyComponentVSG(chrono::vsg3d::ChVisualSystemVSG* vsys, bool& buttonPressed)
        : m_vsys(vsys), pressed(buttonPressed) {}
    virtual void render(vsg::CommandBuffer& cb) override {
        ImGuiWindowFlags window_flags = 0;
        window_flags |= ImGuiWindowFlags_NoTitleBar;
        window_flags |= ImGuiWindowFlags_NoScrollbar;
        window_flags |= ImGuiWindowFlags_NoMove;
        window_flags |= ImGuiWindowFlags_NoResize;
        window_flags |= ImGuiWindowFlags_NoCollapse;
        window_flags |= ImGuiWindowFlags_NoNav;
        window_flags |= ImGuiWindowFlags_NoBackground;

        const ImGuiViewport* viewport = ImGui::GetMainViewport();
        ImGui::SetNextWindowPos(ImVec2(viewport->GetCenter().x, viewport->WorkPos.y + 20), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(300, 0));
        ImGui::Begin("HydroChrono", NULL, window_flags);

        if (ImGui::Button(pressed ? "Playing" : "Paused", ImVec2(200, 40))) {
            pressed = !pressed;
        }

        ImGui::End();
    }

  private:
    chrono::vsg3d::ChVisualSystemVSG* m_vsys;
    bool& pressed;
};

// -----------------------------------------------------------------------------

GUIImplVSG::GUIImplVSG()
    : pVis(chrono_types::make_shared<chrono::vsg3d::ChVisualSystemVSG>()),
      animated_water_(std::make_unique<AnimatedWaterSurface>()) {}

GUIImplVSG::~GUIImplVSG() = default;

void GUIImplVSG::Init(UI& ui, chrono::ChSystem* system, const char* title) {
    system_ = system;  // Cache for time access.

    pVis->AttachSystem(system);

    pVis->SetWindowTitle(title);
    pVis->SetWindowSize(1280, 720);
    pVis->SetWindowPosition(100, 100);
    pVis->SetWindowTitle(title);

    // Deterministic side-on camera: eye on -Y axis, slightly above, looking at origin.
    const chrono::ChVector3d eye(0.0, -kCameraDistance, kCameraHeight);
    const chrono::ChVector3d target(0.0, 0.0, -10.0);
    pVis->AddCamera(eye, target);
    pVis->SetCameraVertical(chrono::CameraVerticalDir::Z);  // enforce Z-up (no roll)
    pVis->SetCameraAngleDeg(40.0);

    // Key light (main directional) — set before Initialize().
    pVis->SetLightIntensity(kKeyIntensity);
    pVis->SetLightDirection(kKeyAzimuth, kKeyElevation);

    // Shadows OFF: Chrono/VSG lacks per-object shadow control. The water plane
    // would cast/receive shadows that flatten craft lighting.
    // pVis->EnableShadows();

    // Skybox provides sky/horizon context, contrasting with dark water plane.
    pVis->EnableSkyBox();

    // Grid disabled: water plane provides the marine context instead.
    // (Keeping this comment for easy re-enable if needed later.)

    if (system) {
        // Water surface will be added when SetWaveModel is called, or a static
        // plane will be added if no wave model is set before first render.
        // For now, defer water surface creation.

        // Apply painted-metal material variants to all existing bodies.
        // Each body gets a subtle color variation for visual distinction.
        int body_index = 0;
        for (auto& body : system->GetBodies()) {
            if (body && body->GetName() != "water_surface" && 
                body->GetName() != "animated_water_surface") {
                const auto material = MakePaintedMetalVariant(body_index);
                ApplyMaterialToAllVisualShapes(*body, material);
                ++body_index;

                // Enable wireframe overlay for improved shape readability (CAD-style).
                if (kEnableWireframe) {
                    auto model = body->GetVisualModel();
                    if (model) {
                        model->EnableWireframe(true);
                    }
                }
            }
        }
    }

    pVis->AddGuiComponent(std::make_shared<MyComponentVSG>(pVis.get(), ui.simulationStarted));

    pVis->Initialize();

    // Add fill light after initialization (scene must exist first).
    AddFillLight(pVis.get());
}

void GUIImplVSG::SetCamera(double x, double y, double z, double dirx, double diry, double dirz) {
    pVis->SetCameraPosition({x, y, z});
    pVis->SetCameraTarget({dirx, diry, dirz});
}

void GUIImplVSG::SetWaveModel(std::shared_ptr<WaveBase> wave) {
    wave_model_ = wave;
    // Water surface will be created/updated on first render via EnsureWaterSurface().
}

void GUIImplVSG::EnsureWaterSurface() {
    // Require both system and visual system to be valid.
    if (!system_ || !pVis) {
        if constexpr (kDebugWaveSurface) {
            static bool printed_once = false;
            if (!printed_once) {
                std::cout << "[WaterSurface] EnsureWaterSurface early return: system_="
                          << (system_ ? "ok" : "null") << " pVis=" << (pVis ? "ok" : "null") << std::endl;
                printed_once = true;
            }
        }
        return;
    }

    // Check if we have a valid wave model (not NoWave).
    bool has_waves = wave_model_ && wave_model_->GetWaveMode() != WaveMode::noWaveCIC;

    // If already initialized for this visual system, nothing to do.
    if (animated_water_->IsInitializedFor(pVis.get())) {
        return;
    }

    if constexpr (kDebugWaveSurface) {
        std::cout << "[WaterSurface] EnsureWaterSurface: has_waves=" << has_waves
                  << " kForceWaterSurface=" << kForceWaterSurface << std::endl;
    }

    // Handle animated or static water surface.
    if (has_waves || kForceWaterSurface) {
        // Initialize animated water surface directly in VSG scene.
        // Must be done after pVis->Initialize() which happens in Init().
        animated_water_->Initialize(pVis.get());

        if (animated_water_->IsInitialized()) {
            // Print status: wave pointer status.
            std::cout << "[WaterSurface] wave_ptr=" << (wave_model_ ? "ok" : "null");
            if (wave_model_) {
                std::cout << " mode=" << static_cast<int>(wave_model_->GetWaveMode());
            }
            std::cout << " " << animated_water_->GetStatusString() << std::endl;

            // Initial update at current system time.
            double t = system_->GetChTime();
            animated_water_->Update(wave_model_, t);
        } else {
            if constexpr (kDebugWaveSurface) {
                std::cout << "[WaterSurface] Initialize() failed!" << std::endl;
            }
        }
    } else if (!water_surface_created_) {
        // No waves and not forcing: add static water plane.
        auto water_plane = CreateStaticWaterPlane();
        system_->AddBody(water_plane);
        water_surface_created_ = true;
        std::cout << "[WaterSurface] static plane created (no wave model)" << std::endl;
    }
}

bool GUIImplVSG::IsRunning(double timestep) {
    if (pVis->Run() == false) return false;

    // Ensure water surface exists (every frame check; returns early if already done).
    EnsureWaterSurface();

    // Update animated water surface if wave model is set (or forced static).
    // VSG vertex buffers are marked dirty() in Update() for GPU re-upload.
    if (animated_water_ && animated_water_->IsInitialized()) {
        double t = system_ ? system_->GetChTime() : 0.0;
        // Update() handles null wave_model_ gracefully (keeps surface flat).
        animated_water_->Update(wave_model_, t);
    }

    pVis->BeginScene();
    pVis->Render();
    pVis->EndScene();

    return true;
}
