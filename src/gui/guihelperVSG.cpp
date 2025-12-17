// =============================================================================
// HydroChrono VSG GUI Implementation
// =============================================================================
// Main implementation of GUIImplVSG - the VSG visualization backend.
// This file coordinates the various VSG subsystems (materials, lighting,
// water surface, GUI overlay) to provide the complete visualization.
// =============================================================================
#include "guihelper_impl.h"

#include "vsg_config.h"
#include "vsg_gui_component.h"
#include "vsg_lighting.h"
#include "vsg_materials.h"
#include "vsg_water_surface.h"

#include <iostream>
#include <memory>

#include <chrono/core/ChTypes.h>

namespace hydroc {
namespace gui {

using namespace vsg_config;

// =============================================================================
// GUIImplVSG Implementation
// =============================================================================

GUIImplVSG::GUIImplVSG()
    : pVis(chrono_types::make_shared<chrono::vsg3d::ChVisualSystemVSG>()),
      animated_water_(std::make_unique<AnimatedWaterSurface>()),
      viewer_settings_(std::make_unique<ViewerSettings>()) {}

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
    pVis->SetCameraVertical(chrono::CameraVerticalDir::Z);  // Enforce Z-up (no roll)
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

    pVis->AddGuiComponent(std::make_shared<HydroChronoGuiComponent>(pVis.get(), ui.simulationStarted,
                                                                    viewer_settings_.get()));

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
        // Use viewer_settings_ resolution if available.
        int resolution = (viewer_settings_) ? viewer_settings_->grid_resolution : 0;
        animated_water_->Initialize(pVis.get(), resolution);

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
    if (pVis->Run() == false) {
        return false;
    }

    // Ensure water surface exists (every frame check; returns early if already done).
    EnsureWaterSurface();

    // Handle viewer settings changes.
    if (viewer_settings_ && animated_water_) {
        // Handle visibility toggle.
        animated_water_->SetVisible(viewer_settings_->show_water);

        // Handle resolution change (requires mesh rebuild).
        if (viewer_settings_->resolution_changed && animated_water_->IsInitialized()) {
            animated_water_->Reinitialize(viewer_settings_->grid_resolution);
            viewer_settings_->resolution_changed = false;
        }
    }

    // Update animated water surface if wave model is set (or forced static).
    // VSG vertex buffers are marked dirty() in Update() for GPU re-upload.
    if (animated_water_ && animated_water_->IsInitialized()) {
        double t = system_ ? system_->GetChTime() : 0.0;
        // Update() handles null wave_model_ gracefully (keeps surface flat).
        // Pass viewer_settings_ for scale multiplier and throttle.
        animated_water_->Update(wave_model_, t, viewer_settings_.get());
    }

    pVis->BeginScene();
    pVis->Render();
    pVis->EndScene();

    return true;
}

}  // namespace gui
}  // namespace hydroc
