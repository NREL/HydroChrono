// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Radu Serban, Zuriah Quinton
// =============================================================================
//
// Recall that Irrlicht uses a left-hand frame, so everything is rendered with
// left and right flipped.
//
// =============================================================================

#include <cstdio>

#include "chrono/solver/ChSolverPMINRES.h"
#include "chrono/assets/ChPointPointDrawing.h"

#include "chrono/physics/ChSystemNSC.h"
#include "chrono/solver/ChIterativeSolverLS.h"
#include "chrono/timestepper/ChTimestepper.h"
#include "chrono/physics/ChBody.h"

#include "chrono_irrlicht/ChIrrApp.h"
#include "chrono_irrlicht/ChIrrMeshTools.h"

#include "chrono/fea/ChElementAddedMass.h"
#include "chrono/fea/ChNodeFEAxyz.h"
#include "chrono/fea/ChLinkPointFrame.h"
#include "chrono/motion_functions/ChFunctionPosition.h"
#include "chrono/physics/ChForce.h"

#include "H5Cpp.h"

using namespace chrono;
using namespace chrono::irrlicht;
using namespace chrono::fea;

using namespace irr;
using namespace irr::core;
using namespace irr::scene;
using namespace irr::video;
using namespace irr::io;
using namespace irr::gui;

// =============================================================================
class MyEventReceiver : public IEventReceiver {
  public:
    MyEventReceiver(ChIrrAppInterface* myapp, bool& buttonPressed)
        : pressed(buttonPressed) {
        // store pointer application
        application = myapp;

        // ..add a GUI button to control pause/play
        pauseButton = application->GetIGUIEnvironment()->addButton(rect<s32>(510, 20, 650, 35));
        buttonText = application->GetIGUIEnvironment()->addStaticText(L"Paused", rect<s32>(560, 20, 600, 35), false);
    }

    bool OnEvent(const SEvent& event) {
        // check if user clicked button
        if (event.EventType == EET_GUI_EVENT) {
            switch (event.GUIEvent.EventType) {
                case EGET_BUTTON_CLICKED:
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

  private:
    ChIrrAppInterface* application;
    IGUIButton* pauseButton;
    IGUIStaticText* buttonText;

    bool& pressed;
};

// =============================================================================
class LinRestorFile {
private:
    ChMatrixDynamic<double> lin_matrix;
    ChVector<double> equil_pos; 
    ChVector<double> equil_rot;
    std::string h5_file_name;
    std::string matrix_data_name;
    std::string cb_data_name;
    std::string cg_data_name;
    bool file_set_up = false;
    bool data_set_up = false;
    bool init_matrix = false;
    void read_data() {
        // testing HDF5 compatibility
        // open file with read only access
        H5::H5File sphereFile(h5_file_name, H5F_ACC_RDONLY); // "../../test_for_chrono/sphere.h5"
        H5::DataSet dataset = sphereFile.openDataSet(matrix_data_name); // "body1/hydro_coeffs/linear_restoring_stiffness"
        // Get filespace for rank and dimension
        H5::DataSpace filespace = dataset.getSpace();
        // Get number of dimensions in the file dataspace
        // Get and print the dimension sizes of the file dataspace
        hsize_t dims[2];    // dataset dimensions
        int rank = filespace.getSimpleExtentDims(dims);
        // read file into data_out 2d array
        H5::DataSpace mspace1(rank, dims);
        double temp[36];
        // read file info into data_out, a 2d array
        dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
        // turn the 2d array into a ChMatrix (Eigen dynamic matrix)
        lin_matrix.resize(6, 6);
        for (int i = 0; i < dims[0]; i++) {
            for (int j = 0; j < dims[1]; j++) {
                lin_matrix(i, j) = temp[i*dims[1]+j];
            }
        }
        dataset.close();

        // repeat same steps from above to get the cb and cg...reusing some of the previous arrays etc
        dataset = sphereFile.openDataSet(cb_data_name);        
        filespace = dataset.getSpace();
        rank = filespace.getSimpleExtentDims(dims);
        mspace1 = H5::DataSpace(rank, dims);
        dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
        // put into equil chvector
        for (int i = 0; i < 3; i++) {
            equil_rot[i] = temp[i];
        }
        dataset.close();

        // repeat finally for cg
        dataset = sphereFile.openDataSet(cg_data_name);
        filespace = dataset.getSpace();
        rank = filespace.getSimpleExtentDims(dims);
        mspace1 = H5::DataSpace(rank, dims);
        dataset.read(temp, H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
        // put into equil chvector
        for (int i = 0; i < 3; i++) {
            equil_pos[i] = temp[i];
        }
        dataset.close();

        sphereFile.close();
        init_matrix = true;
    }

public:
    LinRestorFile() {
        file_set_up = false;
        data_set_up = false;
    }

    LinRestorFile(std::string file, std::string m_data_name, std::string cg, std::string cb) {
        file_set_up = true;
        data_set_up = true;
        h5_file_name = file;
        matrix_data_name = m_data_name;
        cg_data_name = cg;
        cb_data_name = cb;
        read_data();
        //debug statement to check read in variables
        //std::cout << "testing...............\nlinear restoring matrix is\n" << lin_matrix << "\ncg is\n" << equil_pos << "\ncb is\n" << equil_rot << std::endl;
    }
    LinRestorFile(const LinRestorFile& other) {
        file_set_up = other.file_set_up;
        data_set_up = other.data_set_up;
        init_matrix = other.init_matrix;
        if (file_set_up && data_set_up) {
            h5_file_name = other.h5_file_name;
            matrix_data_name = other.matrix_data_name;
        }
        if (init_matrix) {
            lin_matrix = other.lin_matrix;
            equil_pos = other.equil_pos;
        }
    }
    ~LinRestorFile() {}
    void set_file_name(std::string name) {
        h5_file_name = name;
        file_set_up = true;
        if (!init_matrix && file_set_up && data_set_up) {
            read_data();
        }
    }
    void set_data_name(std::string location) {
        matrix_data_name = location;
        data_set_up = true;
        if (!init_matrix && data_set_up && file_set_up) {
            read_data();
        }
    }
    ChMatrixDynamic<double> get_lin_matrix() const {
        return lin_matrix;
    }
    ChVector<> get_equil_cg() const {
        return equil_pos;
    }
    ChVector<> get_equil_cb() const {
        return equil_rot;
    }

};
// =============================================================================
class LinRestorForce {
  private:
    std::shared_ptr<ChBody> bobber;
    LinRestorFile fileInfo;
    ChVectorN<double, 6> equil;

  public: 
    LinRestorForce() {}
    LinRestorForce(LinRestorFile& lin, std::shared_ptr<ChBody> object) {
        bobber = object;
        fileInfo = lin;
        equil << fileInfo.get_equil_cg().eigen(), fileInfo.get_equil_cb().eigen();
    }
    LinRestorForce(const LinRestorForce& other) {
        *this = other;
    } 
    ~LinRestorForce() {} //TODO

    LinRestorForce operator = (const LinRestorForce& rhs) {
        bobber = rhs.bobber;
        fileInfo = rhs.fileInfo;
        equil << fileInfo.get_equil_cg().eigen(), fileInfo.get_equil_cb().eigen();
    }

    ChVectorN<double, 6> Get_p() const {
        ChVectorN<double, 6> temp;
        //ChVector<double> force_or_torque;
        //std::cout << "bobber pos and rot:\n" << bobber->GetPos() << bobber->GetRot().Q_to_Euler123().eigen() << std::endl;
        temp << bobber->GetPos().eigen(), bobber->GetRot().Q_to_Euler123().eigen();
        temp = equil - temp;
        temp = fileInfo.get_lin_matrix() * temp;
        //std::cout << "force should be....\n" << temp << std::endl;
        return temp;
    }

    double coordinateFunc(int i) {
        if (i >= 0 && i < 6) {
            //std::cout << "force still is....\n" << Get_p() << std::endl;
            return Get_p()[i];
        }
        else {
            std::cout << "wrong index" << std::endl;
            return 0;
        }
    }
    //LinRestorForce& set_force_torque(ChForce::ForceType t) {
    //    type = t;
    //    return *this;
    //}
};
// =============================================================================
class ForceTorqueFunc : public ChFunction {
private:
    LinRestorForce& base;
    int index;
public:
    ForceTorqueFunc(LinRestorForce& b, int i) : base(b), index(i) { }
    virtual ForceTorqueFunc* Clone() const override { return new ForceTorqueFunc(*this); }
    virtual double 	Get_y(double x) const override { 
        //std::cout << "getting y..." << std::endl;
        return base.coordinateFunc(index);  }
};

int main(int argc, char* argv[]) {
    GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

    ChSystemNSC system;
    system.Set_G_acc(ChVector<>(0, -9.81, 0));

    // Create the ground body with two visualization spheres
    // -----------------------------------------------------
    double rest_length = 1.5;
    double spring_coef = 50;
    double damping_coef = 0.5;

    auto ground = chrono_types::make_shared<ChBody>();
    system.AddBody(ground);
    ground->SetIdentifier(-1);
    ground->SetBodyFixed(true);
    ground->SetCollide(false);

    {
        auto sph_1 = chrono_types::make_shared<ChSphereShape>();
        sph_1->GetSphereGeometry().rad = 0.1;
        sph_1->Pos = ChVector<>(-1, 0, 0);
        ground->AddAsset(sph_1);

        auto sph_2 = chrono_types::make_shared<ChSphereShape>();
        sph_2->GetSphereGeometry().rad = 0.1;
        sph_2->Pos = ChVector<>(1, 0, 0);
        ground->AddAsset(sph_2);
    }

    // Create a body suspended through a ChLinkTSDA
    // -------------------------------------------------------------

    auto body_1 = chrono_types::make_shared<ChBody>();
    system.AddBody(body_1);
    body_1->SetPos(ChVector<>(-1, -3, 0));
    body_1->SetIdentifier(1);
    body_1->SetBodyFixed(false);
    body_1->SetCollide(false);
    body_1->SetMass(1);
    body_1->SetInertiaXX(ChVector<>(1, 1, 1));

    // Attach a visualization asset.
    auto box_1 = chrono_types::make_shared<ChBoxShape>();
    box_1->GetBoxGeometry().SetLengths(ChVector<>(1, 1, 1));
    body_1->AddAsset(box_1);
    auto col_1 = chrono_types::make_shared<ChColorAsset>();
    col_1->SetColor(ChColor(0.6f, 0, 0));
    body_1->AddAsset(col_1);

    // Create the spring between body_1 and ground. The spring end points are
    // specified in the body relative frames.
    auto spring_1 = chrono_types::make_shared<ChLinkTSDA>();
    spring_1->Initialize(body_1, ground, true, ChVector<>(0, 0, 0), ChVector<>(-1, 0, 0), false, rest_length);
    spring_1->SetSpringCoefficient(spring_coef);
    spring_1->SetDampingCoefficient(damping_coef);
    system.AddLink(spring_1);

    // Attach a visualization asset.
    spring_1->AddAsset(col_1);
    spring_1->AddAsset(chrono_types::make_shared<ChPointPointSpring>(0.05, 80, 15));

    // Create a body suspended through a ChLinkTSDA (now add added mass)
    // -------------------------------------------------------------------

    //auto body_2 = chrono_types::make_shared<ChBody>();
    auto body_2 = chrono_types::make_shared<ChBody>();
    system.AddBody(body_2);
    body_2->SetPos(ChVector<>(1, -3, 0));
    body_2->SetIdentifier(1);
    body_2->SetBodyFixed(false);
    body_2->SetCollide(false);
    body_2->SetMass(1);
    body_2->SetInertiaXX(ChVector<>(1, 1, 1));

    // Attach a visualization asset.
    auto box_2 = chrono_types::make_shared<ChBoxShape>();
    box_2->GetBoxGeometry().SetLengths(ChVector<>(1, 1, 1));
    body_2->AddAsset(box_1);
    auto col_2 = chrono_types::make_shared<ChColorAsset>();
    col_2->SetColor(ChColor(0, 0, 0.6f));
    body_2->AddAsset(col_2);

    auto spring_2 = chrono_types::make_shared<ChLinkTSDA>();
    spring_2->Initialize(body_2, ground, true, ChVector<>(0, 0, 0), ChVector<>(1, 0, 0), false, rest_length);
    spring_2->SetSpringCoefficient(spring_coef);
    spring_2->SetDampingCoefficient(damping_coef);
    system.AddLink(spring_2);

    // Attach a visualization asset.
    spring_2->AddAsset(col_2);
    spring_2->AddAsset(chrono_types::make_shared<ChPointPointSpring>(0.05, 80, 15));

     //now add a dumb node and link it to body_2 for future use with the added_mass element.
     //Node and link must be added to same system as body. Node can only be added through a mesh
    //auto dumbNode = chrono_types::make_shared<ChNodeFEAxyz>(body_2->GetPos());
    //auto mesh = chrono_types::make_shared<ChMesh>();
    //mesh->AddNode(dumbNode);
    //system.Add(mesh);
    //body_2->SetBodyFixed(false);
    //auto nodeLink = chrono_types::make_shared<ChLinkPointFrame>();
    //nodeLink->Initialize(dumbNode, body_2);
    //system.AddLink(nodeLink); 

    //std::vector<std::shared_ptr<ChNodeFEAxyz> > nodesForAddedMass;
    //nodesForAddedMass.push_back(dumbNode);
    //auto added_mass = chrono_types::make_shared<ChElementAddedMass>();
    //added_mass->SetNodes(nodesForAddedMass);
    //mesh->AddElement(added_mass);   

    //ChSparseMatrix M; 
    //system.GetMassMatrix(&M);
    //std::cout << "initial mass matrix\n" << M << std::endl;


    // testing adding external forces to the body
    LinRestorFile test("../../test_for_chrono/sphere.h5", "body1/hydro_coeffs/linear_restoring_stiffness", "body1/properties/cg", "body1/properties/cb");
    LinRestorForce blah(test, body_2);
    auto force = chrono_types::make_shared<ChForce>();
    auto forcex = chrono_types::make_shared<ForceTorqueFunc>(blah, 0);
    force->SetF_x(forcex);
    auto forcey = chrono_types::make_shared<ForceTorqueFunc>(blah, 1);
    force->SetF_y(forcey);
    auto forcez = chrono_types::make_shared<ForceTorqueFunc>(blah, 2);
    force->SetF_z(forcez);
    auto torque = chrono_types::make_shared<ChForce>();
    torque->SetMode(ChForce::ForceType::TORQUE);
    auto torquex = chrono_types::make_shared<ForceTorqueFunc>(blah, 3);
    torque->SetF_x(torquex);
    auto torquey = chrono_types::make_shared<ForceTorqueFunc>(blah, 4);
    torque->SetF_y(torquey);
    auto torquez = chrono_types::make_shared<ForceTorqueFunc>(blah, 5);
    torque->SetF_z(torquez);

    body_2->AddForce(force);

    // Create the Irrlicht application
    // -------------------------------

    ChIrrApp application(&system, L"ChAddedMass Demo", core::dimension2d<u32>(800, 600), VerticalDir::Y);
    application.AddTypicalLogo();
    application.AddTypicalSky();
    application.AddTypicalLights();
    application.AddTypicalCamera(core::vector3df(0, 0, 6));

    application.AssetBindAll();
    application.AssetUpdateAll();

    // some tools to handle the pause button
    bool buttonPressed = false;
    MyEventReceiver receiver(&application, buttonPressed);
    application.SetUserEventReceiver(&receiver);

    
    // Info about which solver to use - may want to change this later
    auto gmres_solver = chrono_types::make_shared<ChSolverMINRES>();  // change to mkl or minres solver? MINRES
                                                                          // requires symmetric matrix
    gmres_solver->SetMaxIterations(300);
    system.SetSolver(gmres_solver);
    double timestep = 0.001;
    application.SetTimestep(timestep);

    /*system.GetMassMatrix(&M);
    std::cout << "initial mass matrix\n" << M << std::endl;*/

    GetLog().SetNumFormat("%10.5f");

    // Simulation loop
    int frame = 0;
    while (application.GetDevice()->run()) {
        application.BeginScene();

        application.DrawAll();
        if (buttonPressed) {
        if (frame == 0) {
        }

            application.DoStep();

            if (frame % 50 == 0) {
                GetLog() << system.GetChTime() << "  " << spring_1->GetLength() << "  " << spring_1->GetVelocity()
                         << "  " << spring_1->GetForce() << "\n";

                GetLog() << "            " << spring_2->GetLength() << "  " << spring_2->GetVelocity() << "  "/*
                         << body_2->GetAppliedForce() << "\n\n"*/;
                /*system.GetMassMatrix(&M);
                std::cout << "mass matrix\n" << M << std::endl;*/
            }
            frame++;
        }

        application.EndScene();
    }

    return 0;
}
