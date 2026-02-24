// =============================================================================
// HydroChrono VSG GUI Component - Implementation
// =============================================================================

// Suppress MSVC warning about IMGUI_API macro redefinition in vsgImGui headers.
// This is a harmless conflict between vsgImGui/imgui.h and vsgImGui/Export.h.
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4005)  // macro redefinition
#endif

#include "vsg_gui_component.h"
#include "vsg_mooring_lines.h"

#include <vsgImGui/imgui.h>

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#include <algorithm>
#include <cstdio>

namespace hydroc {
namespace gui {

HydroChronoGuiComponent::HydroChronoGuiComponent(chrono::vsg3d::ChVisualSystemVSG* vsys, bool& button_pressed,
                                                 ViewerSettings* settings,
                                                 MooringLinesViz* mooring_viz)
    : vsys_(vsys), pressed_(button_pressed), settings_(settings),
      mooring_viz_(mooring_viz) {}

void HydroChronoGuiComponent::render(vsg::CommandBuffer& /*cb*/) {
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

        // ---- Waves ----
        if (ImGui::CollapsingHeader("Waves", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::Checkbox("Show Water", &settings_->show_water);

            ImGui::Checkbox("Show Wireframe", &settings_->show_water_grid);
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Show mesh edges on water surface.\n"
                                  "Helps visualize wave motion and ripples.");
            }

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

            ImGui::Checkbox("Radiation (approx.)", &settings_->show_radiation_viz);
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Visual approximation of radiated waves.\n"
                                  "Does NOT affect physics - for feedback only.");
            }

            if (settings_->show_radiation_viz) {
                ImGui::Indent();
                ImGui::SliderFloat("Visual Scale##rad", &settings_->radiation_visual_scale, 0.1f, 5.0f, "%.1fx");
                if (ImGui::IsItemHovered()) {
                    ImGui::SetTooltip("Visual amplification for inspection.\n"
                                      "1.0x = baseline, increase to see small ripples.\n"
                                      "Does NOT affect physics.");
                }
                ImGui::TextDisabled("Wavelength/speed derived from wave period");
                ImGui::Unindent();
            }
        }

        // ---- Mooring ----
        RenderMooringPanel();

        // ---- Display ----
        if (ImGui::CollapsingHeader("Display")) {
            ImGui::SliderInt("Update Hz", &settings_->update_hz, 10, 60, "%d Hz");

            ImGui::Checkbox("Show Status", &settings_->show_water_status);
            if (settings_->show_water_status) {
                ImGui::TextDisabled("Water: %s | Grid: %d | Hz: %d",
                                    settings_->show_water ? "ON" : "OFF",
                                    settings_->grid_resolution,
                                    settings_->update_hz);
            }
        }

        ImGui::End();
    }

    RenderColorBar();
}

// ---------------------------------------------------------------------------

void HydroChronoGuiComponent::RenderMooringPanel() {
    if (!settings_ || !mooring_viz_ || !mooring_viz_->IsInitialized()) return;

    if (!ImGui::CollapsingHeader("Mooring", ImGuiTreeNodeFlags_DefaultOpen)) return;

    ImGui::Checkbox("Colour by Tension", &settings_->show_mooring_colors);
    if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip("Map tension magnitude to colour\n"
                          "along each mooring line.");
    }

    if (settings_->show_mooring_colors) {
        ImGui::Indent();

        ImGui::Checkbox("Lock Range", &settings_->mooring_range_locked);
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Freeze the colour range for stable comparison.\n"
                              "Uncheck to auto-track.");
        }

        ImGui::TextDisabled("Range: %.3g .. %.3g N",
                            mooring_viz_->ScalarMin(),
                            mooring_viz_->ScalarMax());

        ImGui::Unindent();
    }
}

// ---------------------------------------------------------------------------

void HydroChronoGuiComponent::RenderColorBar() {
    if (!settings_ || !settings_->show_mooring_colors) return;
    if (!mooring_viz_ || !mooring_viz_->IsInitialized()) return;

    const ImGuiViewport* vp = ImGui::GetMainViewport();
    const float bar_w = 20.0f;
    const float bar_h = 200.0f;
    const float margin = 30.0f;
    const float label_gap = 6.0f;

    const float x0 = vp->WorkPos.x + vp->WorkSize.x - margin - bar_w;
    const float y0 = vp->WorkPos.y + vp->WorkSize.y * 0.5f - bar_h * 0.5f;

    ImDrawList* dl = ImGui::GetForegroundDrawList();

    // Draw the gradient as a stack of thin horizontal strips.
    constexpr int kStrips = 64;
    const float strip_h = bar_h / kStrips;
    for (int i = 0; i < kStrips; ++i) {
        const float t_top = 1.0f - static_cast<float>(i) / kStrips;
        const float t_bot = 1.0f - static_cast<float>(i + 1) / kStrips;

        auto to_im = [](float t) -> ImU32 {
            t = std::clamp(t, 0.0f, 1.0f);
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
            return IM_COL32(static_cast<int>(r * 255),
                            static_cast<int>(g * 255),
                            static_cast<int>(b * 255), 255);
        };

        ImU32 col_top = to_im(t_top);
        ImU32 col_bot = to_im(t_bot);
        dl->AddRectFilledMultiColor(
            ImVec2(x0, y0 + i * strip_h),
            ImVec2(x0 + bar_w, y0 + (i + 1) * strip_h),
            col_top, col_top, col_bot, col_bot);
    }

    // Outline.
    dl->AddRect(ImVec2(x0, y0), ImVec2(x0 + bar_w, y0 + bar_h),
                IM_COL32(200, 200, 200, 180));

    // Labels.
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.3g N", mooring_viz_->ScalarMax());
    dl->AddText(ImVec2(x0 - label_gap - ImGui::CalcTextSize(buf).x, y0 - 2.0f),
                IM_COL32(220, 220, 220, 255), buf);

    std::snprintf(buf, sizeof(buf), "%.3g N", mooring_viz_->ScalarMin());
    dl->AddText(ImVec2(x0 - label_gap - ImGui::CalcTextSize(buf).x,
                        y0 + bar_h - ImGui::GetTextLineHeight() + 2.0f),
                IM_COL32(220, 220, 220, 255), buf);

    const char* title = "Tension";
    dl->AddText(ImVec2(x0 + bar_w * 0.5f - ImGui::CalcTextSize(title).x * 0.5f,
                        y0 - ImGui::GetTextLineHeight() - 4.0f),
                IM_COL32(220, 220, 220, 255), title);
}

}  // namespace gui
}  // namespace hydroc

