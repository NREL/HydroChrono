#pragma once

#include <memory>

#include <hydroc/config.h>
#include <hydroc/gui/guihelper.h>
#include <hydroc/logging.h>
#include <hydroc/waves/wave_base.h>

#include <chrono/physics/ChSystem.h>

#include <chrono/core/ChCoordsys.h>
#include <chrono/core/ChQuaternion.h>
#include <chrono/core/ChVector3.h>

#include <chrono/assets/ChVisualSystem.h>

#ifdef HYDROCHRONO_HAVE_IRRLICHT
    #include <chrono_irrlicht/ChVisualSystemIrrlicht.h>
#endif

#ifdef HYDROCHRONO_HAVE_VSG
    #include "chrono_vsg/ChVisualSystemVSG.h"
#endif

namespace chrono {
class ChSystem;
}

namespace hydroc {
namespace gui {

// Base class for private GUI implementation.
class GUIImpl {
  public:
    GUIImpl() {}
    virtual ~GUIImpl() {}
    GUIImpl(const GUIImpl&)            = delete;
    GUIImpl& operator=(const GUIImpl&) = delete;

    virtual void Init(UI& ui, chrono::ChSystem*, const char* title) {
        hydroc::cli::LogWarning(
            "Warning: GUI deactivated. HydroChrono was built without run-time visualization support.");
    }
    virtual void SetCamera(double /*x*/, double /*y*/, double /*z*/,
                           double /*dirx*/, double /*diry*/, double /*dirz*/) {}
    virtual bool IsRunning(double /*timestep*/) { return true; }
    virtual void SetWaveModel(std::shared_ptr<WaveBase> /*wave*/) {}
    virtual void SetWaterGridExtent(double /*width*/, double /*length*/,
                                     double /*center_x*/, double /*center_y*/) {}
    virtual void SetMooringLineProvider(MooringVizProvider /*provider*/) {}
};

#ifdef HYDROCHRONO_HAVE_IRRLICHT
class GUIImplIRR : public GUIImpl {
  public:
    GUIImplIRR();
    GUIImplIRR(const GUIImplIRR&)            = delete;
    GUIImplIRR& operator=(const GUIImplIRR&) = delete;

    virtual void Init(UI& ui, chrono::ChSystem*, const char* title) override;
    virtual void SetCamera(double x, double y, double z, double dirx, double diry, double dirz) override;
    virtual bool IsRunning(double timestep) override;

  private:
    class MyActionReceiver;
    void InitReceiver(bool& simulationStarted);
    std::shared_ptr<chrono::irrlicht::ChVisualSystemIrrlicht> pVis;
    std::shared_ptr<MyActionReceiver> receiver;
};
#endif

#ifdef HYDROCHRONO_HAVE_VSG
// Forward declarations (defined in vsg_water_surface.h, vsg_gui_component.h, vsg_mooring_lines.h).
class AnimatedWaterSurface;
class MooringLinesViz;
struct ViewerSettings;

class GUIImplVSG : public GUIImpl {
  public:
    GUIImplVSG();
    ~GUIImplVSG();
    GUIImplVSG(const GUIImplVSG&)            = delete;
    GUIImplVSG& operator=(const GUIImplVSG&) = delete;

    virtual void Init(UI& ui, chrono::ChSystem*, const char* title) override;
    virtual void SetCamera(double x, double y, double z, double dirx, double diry, double dirz) override;
    virtual bool IsRunning(double timestep) override;
    virtual void SetWaveModel(std::shared_ptr<WaveBase> wave) override;
    virtual void SetWaterGridExtent(double width, double length, double center_x, double center_y) override;
    virtual void SetMooringLineProvider(MooringVizProvider provider) override;

  private:
    /// Ensure water surface is created (animated if wave model set, static otherwise).
    void EnsureWaterSurface();

    /// Update radiation source bodies for visualization.
    /// Iterates over all moving bodies and updates their radiation state.
    void UpdateRadiationSourceBody(double t);

    std::shared_ptr<chrono::vsg3d::ChVisualSystemVSG> pVis;
    std::shared_ptr<WaveBase> wave_model_;
    chrono::ChSystem* system_ = nullptr;
    std::unique_ptr<AnimatedWaterSurface> animated_water_;
    std::unique_ptr<ViewerSettings> viewer_settings_;

    MooringVizProvider mooring_provider_;
    std::unique_ptr<MooringLinesViz> mooring_viz_;
    vsg::ref_ptr<vsg::Group> mooring_scene_group_;
};
#endif

}  // namespace gui
}  // namespace hydroc
