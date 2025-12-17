// =============================================================================
// HydroChrono VSG GUI Component - Implementation
// =============================================================================
#include "vsg_gui_component.h"

#include <vsgImGui/imgui.h>

namespace hydroc {
namespace gui {

HydroChronoGuiComponent::HydroChronoGuiComponent(chrono::vsg3d::ChVisualSystemVSG* vsys, bool& button_pressed,
                                                 ViewerSettings* settings)
    : vsys_(vsys), pressed_(button_pressed), settings_(settings) {}

void HydroChronoGuiComponent::render(vsg::CommandBuffer& cb) {
    // Play/Pause button (top center, no background).
    {
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
        ImGui::Begin("HydroChrono", nullptr, window_flags);

        if (ImGui::Button(pressed_ ? "Playing" : "Paused", ImVec2(200, 40))) {
            pressed_ = !pressed_;
        }

        ImGui::End();
    }

    // Viewer Settings panel (bottom left corner).
    if (settings_) {
        ImGuiWindowFlags settings_flags = 0;
        settings_flags |= ImGuiWindowFlags_AlwaysAutoResize;

        const ImGuiViewport* viewport = ImGui::GetMainViewport();
        ImGui::SetNextWindowPos(ImVec2(viewport->WorkPos.x + 10,
                                       viewport->WorkPos.y + viewport->WorkSize.y - 10),
                                ImGuiCond_FirstUseEver, ImVec2(0.0f, 1.0f));  // Anchor bottom-left

        ImGui::Begin("Viewer Settings", nullptr, settings_flags);

        // Water visibility toggle.
        ImGui::Checkbox("Show Water", &settings_->show_water);

        // Grid resolution combo box.
        {
            const char* resolution_labels[] = {"32x32", "64x64", "96x96", "128x128"};
            int current_idx = 0;
            for (int i = 0; i < ViewerSettings::kResolutionCount; ++i) {
                if (ViewerSettings::kResolutionOptions[i] == settings_->grid_resolution) {
                    current_idx = i;
                    break;
                }
            }
            int prev_idx = current_idx;
            if (ImGui::Combo("Grid Resolution", &current_idx, resolution_labels, ViewerSettings::kResolutionCount)) {
                if (current_idx != prev_idx) {
                    settings_->grid_resolution = ViewerSettings::kResolutionOptions[current_idx];
                    settings_->resolution_changed = true;
                }
            }
        }

        // Update rate slider.
        ImGui::SliderInt("Update Hz", &settings_->update_hz, 10, 60, "%d Hz");

        // Status line toggle.
        ImGui::Checkbox("Show Status", &settings_->show_water_status);

        // Optional status display.
        if (settings_->show_water_status) {
            ImGui::Separator();
            ImGui::TextDisabled("Water: %s | Grid: %d | Hz: %d",
                                settings_->show_water ? "ON" : "OFF",
                                settings_->grid_resolution,
                                settings_->update_hz);
        }

        ImGui::End();
    }
}

}  // namespace gui
}  // namespace hydroc

