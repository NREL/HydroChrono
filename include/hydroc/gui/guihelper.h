#pragma once

#include <memory>

namespace chrono {
class ChSystem;
}

namespace hydroc {
namespace gui {

struct UI {
    UI() {}
    virtual ~UI() {}

    UI(const UI&) = delete;
    UI& operator=(const UI&) = delete;

    /**@brief Initialize the system
     *
     * Should be called after the given ChSystem is fully initialized
     * The best is to call it just before the simulation loop that call IsRunning
     *
     */
    virtual void Init(chrono::ChSystem*, const char* title);

    /**@brief Set Camera position and direction
     *
     */
    virtual void SetCamera(double x, double y, double z, double dirx, double diry, double dirz);

    virtual void ShowFrames(bool show);

    /**@brief To call during simulation loop
     *
     */
    virtual bool IsRunning(double timestep);

    /**@brief return the internal system.
     *
     * Should be called after init.
     */
    chrono::ChSystem* GetSystem() const { return pSystem; }

    bool simulationStarted = true;

  protected:
    chrono::ChSystem* pSystem = nullptr;  // Do not manage the memory
};

class GUIImpl;

struct GUI : public UI {
    GUI();
    GUI(const GUI&) = delete;
    GUI& operator=(const GUI&) = delete;

    void Init(chrono::ChSystem*, const char* title) override;
    void SetCamera(double x, double y, double z, double dirx, double diry, double dirz) override;
    bool IsRunning(double timestep) override;
    void ShowFrames(bool show);

  private:
    std::shared_ptr<hydroc::gui::GUIImpl> pImpl;
};

/**@brief Factory to create UI or GUI
 *
 */
std::shared_ptr<hydroc::gui::UI> CreateUI(bool visualizationOn = true);

}  // end namespace gui
}  // end namespace hydroc