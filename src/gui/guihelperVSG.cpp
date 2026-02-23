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
#include "vsg_mooring_lines.h"
#include "vsg_radiation_surface.h"
#include "vsg_water_surface.h"

#include <hydroc/waves/regular_wave.h>
#include <hydroc/waves/irregular_wave.h>

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
    pVis->EnableSkyTexture(chrono::SkyMode::BOX);

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

void GUIImplVSG::SetWaterGridExtent(double width, double length, double center_x, double center_y) {
    if (!viewer_settings_) {
        return;
    }

    viewer_settings_->grid_width = width;
    viewer_settings_->grid_length = length;
    viewer_settings_->grid_center_x = center_x;
    viewer_settings_->grid_center_y = center_y;
    viewer_settings_->grid_extent_changed = true;

    std::cout << "[WaterSurface] Grid extent: " << width << " x " << length << " m"
              << ", center: (" << center_x << ", " << center_y << ")" << std::endl;
}

void GUIImplVSG::SetMooringLineProvider(MooringVizProvider provider) {
    mooring_provider_ = std::move(provider);
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

    // If already initialized for this visual system, nothing to do.
    if (animated_water_->IsInitializedFor(pVis.get())) {
        return;
    }

    if constexpr (kDebugWaveSurface) {
        bool has_waves = wave_model_ && wave_model_->GetWaveMode() != WaveMode::noWaveCIC;
        std::cout << "[WaterSurface] EnsureWaterSurface: has_waves=" << has_waves << std::endl;
    }

    // Always create animated water surface (supports waves + radiation viz).
    // Static plane fallback removed - animated water handles all cases.
    int resolution = (viewer_settings_) ? viewer_settings_->grid_resolution : 0;
    animated_water_->Initialize(pVis.get(), resolution, viewer_settings_.get());

    if (animated_water_->IsInitialized()) {
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
}

bool GUIImplVSG::IsRunning(double timestep) {
    if (pVis->Run() == false) {
        return false;
    }

    // Reserve a scene-graph group for mooring geometry before the water
    // surface is created, so opaque mooring lines are rendered first and
    // remain visible through the transparent free surface.
    if (mooring_provider_ && !mooring_scene_group_) {
        if (auto scene = pVis->GetVSGScene()) {
            mooring_scene_group_ = vsg::Group::create();
            scene->addChild(mooring_scene_group_);
        }
    }

    EnsureWaterSurface();

    // Mooring line visualisation (lazy-initialised on first valid data).
    if (mooring_provider_) {
        auto line_data = mooring_provider_();
        if (!line_data.empty()) {
            if (!mooring_viz_)
                mooring_viz_ = std::make_unique<MooringLinesViz>();
            if (!mooring_viz_->IsInitializedFor(pVis.get()))
                mooring_viz_->Initialize(pVis.get(), line_data, mooring_scene_group_);
            mooring_viz_->Update(line_data);
        }
    }

    // Handle viewer settings changes.
    if (viewer_settings_ && animated_water_) {
        // Handle visibility toggle.
        animated_water_->SetVisible(viewer_settings_->show_water);

        // Handle resolution or extent change (requires mesh rebuild).
        if ((viewer_settings_->resolution_changed || viewer_settings_->grid_extent_changed) &&
            animated_water_->IsInitialized()) {
            animated_water_->Reinitialize(viewer_settings_->grid_resolution, viewer_settings_.get());
            viewer_settings_->resolution_changed = false;
            viewer_settings_->grid_extent_changed = false;
        }
    }

    // Update animated water surface if wave model is set (or forced static).
    // VSG vertex buffers are marked dirty() in Update() for GPU re-upload.
    if (animated_water_ && animated_water_->IsInitialized()) {
        double t = system_ ? system_->GetChTime() : 0.0;

        // Update radiation source body state if radiation viz is enabled.
        // Chooses the first non-water body in the system as the source.
        // NOTE: This is visualization-only and does NOT affect physics.
        if (viewer_settings_ && viewer_settings_->show_radiation_viz && system_) {
            UpdateRadiationSourceBody(t);
        }

        // Update() handles null wave_model_ gracefully (keeps surface flat).
        // Pass viewer_settings_ for scale multiplier and throttle.
        animated_water_->Update(wave_model_, t, viewer_settings_.get());
    }

    pVis->BeginScene();
    pVis->Render();
    pVis->EndScene();

    return true;
}

void GUIImplVSG::UpdateRadiationSourceBody(double t) {
    if (!animated_water_ || !system_) {
        return;
    }

    // Update global params.
    if (viewer_settings_) {
        RadiationSurfaceViz::Params rad_params;
        rad_params.visual_scale = static_cast<double>(viewer_settings_->radiation_visual_scale);
        rad_params.wave_period = 8.0;  // Default; should match body oscillation frequency

        // Get wave properties from wave model if available.
        if (wave_model_) {
            // Water depth and gravity (available in all wave models).
            if (wave_model_->GetWaterDepth() > 0.0) {
                rad_params.water_depth = wave_model_->GetWaterDepth();
            }
            rad_params.gravity = wave_model_->GetGravity();

            // Wave period from RegularWave.
            if (wave_model_->GetWaveMode() == WaveMode::regular) {
                auto* reg_wave = dynamic_cast<RegularWave*>(wave_model_.get());
                if (reg_wave && reg_wave->regular_wave_omega_ > 0.0) {
                    rad_params.wave_period = (2.0 * M_PI) / reg_wave->regular_wave_omega_;
                }
            }
            // Peak period from IrregularWaves (JONSWAP).
            else if (wave_model_->GetWaveMode() == WaveMode::irregular) {
                auto* irreg_wave = dynamic_cast<IrregularWaves*>(wave_model_.get());
                if (irreg_wave) {
                    // Note: IrregularWaves doesn't expose params_ directly.
                    // For now, keep default. Could add a GetPeakPeriod() accessor.
                }
            }
            // NoWave mode (decay tests): wave_period should ideally be set to the 
            // body's natural period T_n ≈ 2π√((M+A_∞)/K_hs). Using default 8s as fallback.
        }

        animated_water_->GetRadiationViz().SetParams(rad_params);
    }

    // Log source bodies once.
    static bool logged_sources = false;

    // Iterate over ALL moving bodies and update their radiation state.
    for (auto& body : system_->GetBodies()) {
        if (!body) {
            continue;
        }

        const std::string& name = body->GetName();

        // Skip water surfaces and ground bodies.
        if (name == "water_surface" || name == "animated_water_surface" ||
            name == "ground" || name == "floor" || name.empty()) {
            continue;
        }

        // Skip fixed bodies.
        if (body->IsFixed()) {
            continue;
        }

        // Get body motion state.
        const chrono::ChVector3d pos = body->GetPos();
        const chrono::ChVector3d vel = body->GetPosDt();
        const chrono::ChVector3d ang_vel = body->GetAngVelLocal();

        // Estimate body radius from AABB (rough approximation).
        double radius = 5.0;  // Default
        auto aabb = body->GetTotalAABB();
        if (aabb.max.x() > aabb.min.x()) {
            double dx = aabb.max.x() - aabb.min.x();
            double dy = aabb.max.y() - aabb.min.y();
            radius = std::max(dx, dy) / 2.0;
            radius = std::max(radius, 1.0);  // Minimum 1m
        }

        // Update radiation viz for this body.
        animated_water_->GetRadiationViz().SetSourceState(name, pos, vel, ang_vel, radius, t);

        if (!logged_sources) {
            std::cout << "[RadiationViz] Source: " << name << " (r=" << radius << "m)" << std::endl;
        }
    }

    logged_sources = true;
}

}  // namespace gui
}  // namespace hydroc
