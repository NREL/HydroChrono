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
    virtual void SetCamera(double x, double y, double z, double dirx, double diry, double dirz) {}
    virtual bool IsRunning(double timestep) { return true; }
    virtual void SetWaveModel(std::shared_ptr<WaveBase> wave) { (void)wave; }
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
// Forward declaration for animated water surface (defined in guihelperVSG.cpp).
class AnimatedWaterSurface;

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

  private:
    /// Ensure water surface is created (animated if wave model set, static otherwise).
    void EnsureWaterSurface();

    std::shared_ptr<chrono::vsg3d::ChVisualSystemVSG> pVis;
    std::shared_ptr<WaveBase> wave_model_;               ///< Wave model for animated surface
    chrono::ChSystem* system_ = nullptr;                 ///< Cached system pointer for time access
    std::unique_ptr<AnimatedWaterSurface> animated_water_;  ///< Animated water surface (owned)
    bool water_surface_created_ = false;                 ///< True if any water surface was added
};
#endif

}  // namespace gui
}  // namespace hydroc
