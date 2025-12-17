#include "guihelper_impl.h"

#include <algorithm>

#include <chrono/physics/ChBodyEasy.h>

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
// Dark ocean blue-green — darker than sky for clear horizon.
constexpr float kWaterR = 0.03f;
constexpr float kWaterG = 0.10f;
constexpr float kWaterB = 0.14f;

// Semi-transparent surface (60% opacity for visible waterline).
constexpr float kWaterOpacity = 0.60f;

// Subtle specular highlight for calm water.
constexpr float kWaterSpecular = 0.4f;
constexpr float kWaterSpecularExp = 96.0f;

// Slightly glossy, mostly dielectric.
constexpr float kWaterRoughness = 0.15f;
constexpr float kWaterMetallic = 0.0f;

// Water plane geometry.
constexpr double kWaterPlaneSize = 200.0;       // meters (large enough to fill most views)
constexpr double kWaterPlaneThickness = 0.02;   // meters (thin slab)
constexpr double kWaterPlaneZ = -0.05;          // below z=0 to avoid z-fighting with craft

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

chrono::ChVisualMaterial MakeWaterSurfaceMaterial() {
    chrono::ChVisualMaterial mat;

    mat.SetDiffuseColor(chrono::ChColor(kWaterR, kWaterG, kWaterB));
    mat.SetSpecularColor(chrono::ChColor(kWaterSpecular, kWaterSpecular, kWaterSpecular));
    mat.SetSpecularExponent(kWaterSpecularExp);
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
std::shared_ptr<chrono::ChBody> CreateWaterPlane() {
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

GUIImplVSG::GUIImplVSG() : pVis(chrono_types::make_shared<chrono::vsg3d::ChVisualSystemVSG>()) {}

void GUIImplVSG::Init(UI& ui, chrono::ChSystem* system, const char* title) {
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
        // Add water surface plane (visual-only, fixed, non-colliding).
        auto water_plane = CreateWaterPlane();
        system->AddBody(water_plane);

        // Apply painted-metal material variants to all existing bodies.
        // Each body gets a subtle color variation for visual distinction.
        int body_index = 0;
        for (auto& body : system->GetBodies()) {
            if (body && body->GetName() != "water_surface") {
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

bool GUIImplVSG::IsRunning(double timestep) {
    if (pVis->Run() == false) return false;

    pVis->BeginScene();
    pVis->Render();
    pVis->EndScene();

    return true;
}
