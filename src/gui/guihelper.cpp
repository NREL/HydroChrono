#include <hydroc/gui/guihelper.h>
using namespace hydroc::gui;

#include <chrono/physics/ChSystem.h>

#ifdef HYDROCHRONO_HAVE_IRRLICHT
    #include <IEventReceiver.h>  // irrlicht

    #include <chrono_irrlicht/ChIrrMeshTools.h>
    #include <chrono_irrlicht/ChVisualSystemIrrlicht.h>

    #include <chrono/core/ChCoordsys.h>
    #include <chrono/core/ChMathematics.h>
    #include <chrono/core/ChQuaternion.h>
    #include <chrono/core/ChVector.h>

    #include <chrono/assets/ChVisualSystem.h>
    #include <chrono_irrlicht/ChIrrMeshTools.h>
    #include <chrono_irrlicht/ChVisualSystemIrrlicht.h>
#endif

void UI::Init(chrono::ChSystem* system, const char* title) {
    pSystem = system;
}

void UI::SetCamera(double x, double y, double z, double dirx, double diry, double dirz) {}

bool UI::IsRunning(double timestep) {
    return true;
}

std::shared_ptr<hydroc::gui::UI> hydroc::gui::CreateUI(bool visualizationOn) {
    if (visualizationOn) {
        return std::make_shared<hydroc::gui::GUI>();
    } else {
        return std::make_shared<hydroc::gui::UI>();
    }
}

// Private implementation
class hydroc::gui::GUIImpl {
  public:
    GUIImpl();
    GUIImpl(const GUIImpl&)            = delete;
    GUIImpl& operator=(const GUIImpl&) = delete;

    void Init(UI& ui, chrono::ChSystem*, const char* title);
    void SetCamera(double x, double y, double z, double dirx, double diry, double dirz);
    bool IsRunning(double timestep);

  private:
#ifdef HYDROCHRONO_HAVE_IRRLICHT
    class MyActionReceiver;
    void InitReceiver(bool& simulationStarted);
    std::shared_ptr<chrono::irrlicht::ChVisualSystemIrrlicht> pVis;
    std::shared_ptr<MyActionReceiver> receiver;
#endif
};

#ifdef HYDROCHRONO_HAVE_IRRLICHT

using namespace chrono::irrlicht;

using irr::EEVENT_TYPE;
using irr::s32;
using irr::core::rect;
using irr::gui::EGUI_EVENT_TYPE;

///@brief Define a class to manage user inputs via the GUI (i.e. play/pause button)
class hydroc::gui::GUIImpl::MyActionReceiver : public irr::IEventReceiver {
  public:
    MyActionReceiver(bool& buttonPressed);
    bool OnEvent(const irr::SEvent& event);
    void Init(chrono::irrlicht::ChVisualSystemIrrlicht* vsys);
    void SetCamera(double x, double y, double z, double dirx, double diry, double dirz);

  private:
    chrono::irrlicht::ChVisualSystemIrrlicht* vis;
    irr::gui::IGUIButton* pauseButton;
    irr::gui::IGUIStaticText* buttonText;

    bool& pressed;
};

GUIImpl::GUIImpl() : pVis(chrono_types::make_shared<chrono::irrlicht::ChVisualSystemIrrlicht>()) {}

void GUIImpl::InitReceiver(bool& theSimulationStarted) {
    receiver = std::make_shared<MyActionReceiver>(theSimulationStarted);
}

void GUIImpl::Init(UI& ui, chrono::ChSystem* system, const char* title) {
    pVis->AttachSystem(system);

    pVis->SetWindowSize(1280, 720);
    pVis->SetWindowTitle(title);
    pVis->SetCameraVertical(chrono::CameraVerticalDir::Z);
    pVis->Initialize();

    InitReceiver(ui.simulationStarted);
    receiver->Init(pVis.get());
    pVis->AddUserEventReceiver(receiver.get());

    pVis->AddLogo();
    pVis->AddSkyBox();
    pVis->AddCamera(chrono::ChVector<>(8, -25, 15), chrono::ChVector<>(0, 0, 0));
    pVis->AddTypicalLights();
}

void GUIImpl::SetCamera(double x, double y, double z, double dirx, double diry, double dirz) {
    pVis->AddCamera({x, y, z}, {dirx, diry, dirz});
}

bool GUIImpl::IsRunning(double timestep) {
    if (pVis->Run() == false) return false;

    pVis->BeginScene();
    pVis->Render();

    // Add grid to materialize horizontal plane
    tools::drawGrid(pVis.get(), 1, 1, 30, 30,
                    chrono::ChCoordsys<>(chrono::ChVector<>(0, 0.0, 0), chrono::Q_from_AngZ(chrono::CH_C_PI_2)),
                    chrono::ChColor(.1f, .1f, .1f), true);

    pVis->EndScene();
    return true;
}

hydroc::gui::GUIImpl::MyActionReceiver::MyActionReceiver(bool& buttonPressed) : pressed(buttonPressed) {}

/// @brief Initialize Action with System
/// @param vsys
void hydroc::gui::GUIImpl::MyActionReceiver::Init(chrono::irrlicht::ChVisualSystemIrrlicht* vsys) {
    // store pointer application
    vis = vsys;

    // ..add a GUI button to control pause/play
    pauseButton = vis->GetGUIEnvironment()->addButton(rect<s32>(510, 20, 650, 35));
    buttonText  = vis->GetGUIEnvironment()->addStaticText(L"Paused", rect<s32>(560, 20, 600, 35), false);
}

bool hydroc::gui::GUIImpl::MyActionReceiver::OnEvent(const irr::SEvent& event) {
    // check if user clicked button
    if (event.EventType == EEVENT_TYPE::EET_GUI_EVENT) {
        switch (event.GUIEvent.EventType) {
            case EGUI_EVENT_TYPE::EGET_BUTTON_CLICKED:
                pressed = !pressed;
                if (pressed) {
                    buttonText->setText(L"Playing");
                } else {
                    buttonText->setText(L"Paused");
                }
                return pressed;
                break;
            default:
                break;
        }
    }
    return false;
}

#else  // HYDROCHRONO_HAVE_IRRLICHT

GUIImpl::GUIImpl() {}

void GUIImpl::Init(UI& ui, chrono::ChSystem* system, const char* title) {
    std::cout << "Warning: GUI deactivated. Compilation without Irrlicht library" << std::endl;
}

void GUIImpl::SetCamera(double x, double y, double z, double dirx, double diry, double dirz) {}

bool GUIImpl::IsRunning(double timestep) {
    return true;
}

#endif  // HYDROCHRONO_HAVE_IRRLICHT

//

GUI::GUI() : pImpl(std::make_shared<hydroc::gui::GUIImpl>()) {
    simulationStarted = false;  // Simulation is Paused
}

void GUI::Init(chrono::ChSystem* system, const char* title) {
    UI::Init(system, title);
    pImpl->Init(*this, system, title);
}

void GUI::SetCamera(double x, double y, double z, double dirx, double diry, double dirz) {
    pImpl->SetCamera(x, y, z, dirx, diry, dirz);
}

bool GUI::IsRunning(double timestep) {
    if (pImpl->IsRunning(timestep) == false) return false;

    // If still running call the base class
    return UI::IsRunning(timestep);
}
