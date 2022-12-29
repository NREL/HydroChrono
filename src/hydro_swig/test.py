import sys
import math
import pychrono as chrono
import pychrono.irrlicht as chronoirr
import _hydro_py

# create the system
system = chrono.ChSystemNSC()
system.Set_G_acc(chrono.ChVectorD(0, 0, 0))

# use _hydro_py module global function and enum
_hydro_py.test()
b = _hydro_py.NONE
if(b == _hydro_py.NONE):
    print("enum worked")

input = _hydro_py.HydroInputs()
input.regular_wave_amplitude = 1
input.regular_wave_omega = 0.5
print("wave amp = ", input.regular_wave_amplitude)
print("wave omega = ", input.regular_wave_omega)
input.test2()

# define the sphere rigid body from mesh
body = chrono.ChBodyEasyMesh("../../demos/sphere/geometry/oes_task10_sphere.obj",1000,False,True,False)
body.SetMass(261.8e3)
body.SetPos(chrono.ChVectorD(0, 0, -1))
body.SetNameString("body1")
system.AddBody(body)
bodies = _hydro_py.vector_ChBody_ptr(1)
print(bodies.size())
temp = _hydro_py.upcast(body)
# bodies = [body]
# bodies = (body,) # make a list of all bodies (just one here but list all with commas-for c++ maybe look at span)
_hydro_py.TestHydro(_hydro_py.upcast(bodies), "../../demos/sphere/hydroData/sphere.h5", input)

# visualization
irr_app = chronoirr.ChVisualSystemIrrlicht()
irr_app.AttachSystem(system)
irr_app.SetWindowSize(1280, 720)
irr_app.SetWindowTitle("Sphere Decay - hydro_py")
irr_app.SetCameraVertical(chrono.CameraVerticalDir_Z)
irr_app.Initialize()
# irr_app.AddLogo()
irr_app.AddSkyBox("../../data/skybox/")
irr_app.AddCamera(chrono.ChVectorD(0,-30,0), chrono.ChVectorD(0,0,0))
irr_app.AddTypicalLights()


with open("py_out.txt", mode='w', encoding='utf-8') as f:
    # simulation loop
    f.write("#Time	Body Pos	Body vel (heave)	force (heave)")
    # irr_app.SetTimestep(0.001)
    while (irr_app.Run and system.GetChTime() <= 40) :
        irr_app.BeginScene()
        irr_app.Render()
        irr_app.EndScene()
        f.write("\n")
        f.write(str(system.GetChTime()))
        f.write("\t")
        f.write(str(body.GetPos().z))
        system.DoStepDynamics(0.015)
