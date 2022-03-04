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

#include "H5_force_classes.h"
#include "chrono_irrlicht/ChIrrNodeAsset.h"
#include "chrono/physics/ChLoadContainer.h"
#include "chrono/physics/ChLoadsBody.h"

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
				}
				else {
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
class ChLoadAddedMass : public ChLoadCustom {
public:
	ChLoadAddedMass(std::shared_ptr<ChBody> body,  ///< object to apply additional inertia to
		const ChVector<>& m_offset,     ///< offset of the center of mass, in body coordinate system
		const double m_mass/*,            ///< added mass [kg]
		const ChVector<>& m_IXX = VNULL,  ///< added diag. inertia values Ixx, Iyy, Izz (in body coordinate system, centered in body)
		const ChVector<>& m_IXY = VNULL */  ///< added off.diag. inertia values Ixy, Ixz, Iyz including the "-"sign (in body coordinate system, centered in body)
	);

	/// "Virtual" copy constructor (covariant return type).
	virtual ChLoadAddedMass* Clone() const override { return new ChLoadAddedMass(*this); }

	/// Set the inertia tensor of the body, assumed in the body reference basis, with body reference as center.
	/// The provided 3x3 matrix should be symmetric and contain the inertia tensor as:
	/// system: 
	/// <pre>
	///               [ int{y^2+z^2}dm    -int{xy}dm    -int{xz}dm    ]
	/// newXInertia = [                  int{x^2+z^2}   -int{yz}dm    ]
	///               [     (symm.)                    int{x^2+y^2}dm ]
	/// </pre>
	//void SetInertia(const ChMatrix33<>& newXInertia);

	/// Set the inertia tensor of the body, assumed in the body reference basis, with body reference as center.
	/// The return 3x3 symmetric matrix contains the following values:
	/// <pre>
	///  [ int{y^2+z^2}dm    -int{xy}dm    -int{xz}dm    ]
	///  [                  int{x^2+z^2}   -int{yz}dm    ]
	///  [       (symm.)                  int{x^2+y^2}dm ]
	/// </pre>
	//const ChMatrix33<>& GetInertia() const { return this->I; }

	/// Set the diagonal part of the inertia tensor (Ixx, Iyy, Izz values). 
	/// The vector should contain these moments of inertia, assumed in the body reference basis, with body reference as center:
	/// <pre>
	/// iner = [  int{y^2+z^2}dm   int{x^2+z^2}   int{x^2+y^2}dm ]
	/// </pre>
	//void SetInertiaXX(const ChVector<>& iner);

	/// Get the diagonal part of the inertia tensor (Ixx, Iyy, Izz values). 
	/// The vector contains these values, assumed in the body reference basis, with body reference as center:
	/// <pre>
	/// [  int{y^2+z^2}dm   int{x^2+z^2}   int{x^2+y^2}dm ]
	/// </pre>
	//ChVector<> GetInertiaXX() const;

	/// Set the off-diagonal part of the inertia tensor (Ixy, Ixz, Iyz values).
	/// The vector contains these values, assumed in the body reference basis, with body reference as center:
	/// <pre>
	/// iner = [ -int{xy}dm   -int{xz}dm   -int{yz}dm ]
	/// </pre>
	//void SetInertiaXY(const ChVector<>& iner);

	/// Get the extra-diagonal part of the inertia tensor (Ixy, Ixz, Iyz values).
	/// The vector contains these values, assumed in the body reference basis, with body reference as center:
	/// <pre>
	/// [ -int{xy}dm   -int{xz}dm   -int{yz}dm ]
	/// </pre>
	//ChVector<> GetInertiaXY() const;

	/// Compute Q, the generalized load. 
	/// In this case, it computes the quadratic (centrifugal, gyroscopic) terms. 
	/// Signs are negative as Q assumed at right hand side, so Q= -Fgyro -Fcentrifugal
	/// Called automatically at each Update().
	/// The M*a term is not added: to this end one could use LoadIntLoadResidual_Mv afterward.
	virtual void ComputeQ(ChState* state_x,      ///< state position to evaluate Q
		ChStateDelta* state_w  ///< state speed to evaluate Q
	) override;

	/// For efficiency reasons, do not let the parent class do automatic differentiation
	/// to compute the R, K matrices. Use analytic expressions instead. For example, R is 
	/// the well known gyroscopic damping matrix. Also, compute the M matrix.
	virtual void ComputeJacobian(ChState* state_x,       ///< state position to evaluate jacobians
		ChStateDelta* state_w,  ///< state speed to evaluate jacobians
		ChMatrixRef mK,         ///< result -dQ/dx
		ChMatrixRef mR,         ///< result -dQ/dv
		ChMatrixRef mM          ///< result -dQ/da
	) override;

	/// Just for efficiency, override the default LoadIntLoadResidual_Mv, because we can do this in a simplified way.
	virtual void LoadIntLoadResidual_Mv(ChVectorDynamic<>& R,           ///< result: the R residual, R += c*M*w
		const ChVectorDynamic<>& w,     ///< the w vector
		const double c) override;       ///< a scaling factor
private:
	ChVector<> c_m;       ///< offset of center of mass
	double  mass;         ///< added mass
	ChMatrix33<> I;       ///< added inertia tensor, in body coordinates

	virtual bool IsStiff() override { return true; } // this to force the use of the inertial M, R and K matrices

	//static bool use_inertial_damping_matrix_R;  // default true. Can be disabled globally, for testing or optimization
	//static bool use_inertial_stiffness_matrix_K;// default true. Can be disabled globally, for testing or optimization
	//static bool use_gyroscopic_torque;          // default true. Can be disabled globally, for testing or optimization
};
// =============================================================================
int main(int argc, char* argv[]) {
	GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

	ChSystemNSC system;
	system.Set_G_acc(ChVector<>(0, 0, -9.81));

	// Create the Irrlicht application for visualizing
	// names of some functions changed feb 7, 2022 (AddTypicalLogo -> AddLogo, etc)
	ChIrrApp application(&system, L"Sphere Decay Test", core::dimension2d<u32>(800, 600), VerticalDir::Z);
	application.AddLogo();
	application.AddSkyBox();
	application.AddTypicalLights();
	application.AddCamera(core::vector3df(0, 30, 0), core::vector3df(0, 0, 0)); // arguments are (location, orientation) as vectors

	// set up body from a mesh
	std::shared_ptr<ChBody> body = chrono_types::make_shared<ChBodyEasyMesh>(                   //
		GetChronoDataFile("../../test_for_chrono/oes_task10_sphere.obj").c_str(),                 // file name
		1000,                                                                                     // density
		false,                                                                                    // do not evaluate mass automatically
		true,                                                                                     // create visualization asset
		false,                                                                                    // do not collide
		nullptr,                                                                                  // no need for contact material
		0                                                                                         // swept sphere radius
		);

	// old shpere stuff (not mesh)
	//std::shared_ptr<ChBody> body = chrono_types::make_shared<ChBodyEasySphere>(5, 1);
	//auto sph = chrono_types::make_shared<ChSphereShape>();
	//body->AddAsset(sph);

	system.Add(body);
	body->SetPos(ChVector<>(0, 0, -1));
	body->SetMass(261.8e3);
	//body->SetMass(261.8e3 + 130.8768);
	//body->SetMass(261.8e3 + 130.8768e3); // added mass 130.8768 kg times rho=1000
	// attach color asset to body
	auto col_2 = chrono_types::make_shared<ChColorAsset>();
	col_2->SetColor(ChColor(0, 0, 0.6f));
	body->AddAsset(col_2);

	// test added mass as load here-------------------------------------------------------
	//auto my_loadcontainer = chrono_types::make_shared<ChLoadContainer>();

	//vars for testing
	ChVector<> m_offset(0, 0, 0); ///< offset of the center of mass, in body coordinate system
	double m_mass = 130.8768e3; ///< added mass [kg]
	// ChVector<> m_IXX = VNULL;  ///< added diag. inertia values Ixx, Iyy, Izz (in body coordinate system, centered in body)
	// ChVector<> m_IXY = VNULL;   ///< added off.diag. inertia values Ixy, Ixz, Iyz including the "-"sign (in body coordinate system, centered in body)

	auto my_loadcontainer = chrono_types::make_shared< ChLoadContainer>();
	system.Add(my_loadcontainer);
	auto my_loadbodyinertia = chrono_types::make_shared<ChLoadAddedMass>(body, m_offset, m_mass);
	my_loadcontainer->Add(my_loadbodyinertia);

	// testing adding other external forces to the body------------------------------------
	BodyFileInfo sphere_file_info("../../test_for_chrono/sphere.h5", "body1");
	LinRestorForce lin_restor_force_2(sphere_file_info, body);
	ImpulseResponseForce irf(sphere_file_info, body);
	// declare some forces to be initialized in lin_restor_force_2 to be applied to to body later
	auto force = chrono_types::make_shared<ChForce>();
	auto torque = chrono_types::make_shared<ChForce>();
	auto force2 = chrono_types::make_shared<ChForce>();
	auto torque2 = chrono_types::make_shared<ChForce>();
	// set torque flag for torque
	torque->SetMode(ChForce::ForceType::TORQUE);
	torque2->SetMode(ChForce::ForceType::TORQUE);
	// initialize force and torque with member functions
	lin_restor_force_2.SetForce(force);
	lin_restor_force_2.SetTorque(torque);
	irf.SetForce(force2);
	irf.SetTorque(torque2);
	// apply force and torque to the body
	body->AddForce(force);
	body->AddForce(torque);
	body->AddForce(force2);
	body->AddForce(torque2);

	// add buoyancy force from h5 file info
	auto fb = chrono_types::make_shared<BuoyancyForce>(sphere_file_info);
	body->AddForce(fb->getForce_ptr());

	// update irrlicht app with body info
	//application.AssetBind(body_1);
	application.AssetBindAll();
	application.AssetUpdateAll();

	// some tools to handle the pause button
	bool buttonPressed = false;
	MyEventReceiver receiver(&application, buttonPressed);
	application.SetUserEventReceiver(&receiver);

	// Info about which solver to use - may want to change this later
	auto gmres_solver = chrono_types::make_shared<ChSolverMINRES>();  // change to mkl or minres?
	gmres_solver->SetMaxIterations(300);
	system.SetSolver(gmres_solver);
	double timestep = 0.015; // also sets the timesteps in system it seems
	application.SetTimestep(timestep);

	// set up output file for body position each step
	std::ofstream zpos("outfile/addedMassOut.txt", std::ofstream::out);
	zpos.precision(10);
	zpos.width(12);
	//zpos.SetNumFormat("%10.5f"); // 10 characters displayed total, 5 digit precision (after decimal)

	// Simulation loop
	int frame = 0;
	//bool full_period = false;
	//ChVector<> initial_pos = body->GetPos();
	std::cout << "Body mass=" << body->GetMass() << std::endl;
	zpos << "#Time\tBody Pos\tBody vel (heave)\tforce (heave)\n";
	//system.EnableSolverMatrixWrite(true, "blah2");
	while (application.GetDevice()->run() && system.GetChTime() <= 25) {
		application.BeginScene();
		application.DrawAll();
		if (buttonPressed)/*if(true)*/ {
			/*if (frame == 8) {
				ChSparseMatrix M;
				system.GetMassMatrix(&M);
				std::cout << "initial mass matrix\n" << M << std::endl;
				system.DumpSystemMatrices(true, true, true, true, "C:\\Users\\ZQUINTON\\code\\test_for_chrono_build\\Release\\outfile\\addedM");
			}*/
			zpos << system.GetChTime() << "\t" << body->GetPos().z() << "\t" << body->GetPos_dt().z() << "\t"  << body->GetAppliedForce().z() << "\n";
			application.DoStep();
			frame++;
			//if (!full_period && (body->GetPos().Equals(initial_pos, 0.00001) ) && frame > 5 ) {
			//	full_period = true;
			//	std::cout << "frame: " << frame << std::endl;
			//}
			//GetLog() << "\n" << fb->getForce_ptr()->GetForce() << "that's the force pointer value\n\n";
		}

		application.EndScene();
	}
	zpos.close();
	return 0;
}


// -----------------------------------------------------------------------------
// ChLoadAddedMass
// -----------------------------------------------------------------------------

// for testing and optimizations
//bool ChLoadAddedMass::use_inertial_damping_matrix_R = true;   // default true. Can be disabled globally, for testing or optimization
//bool ChLoadAddedMass::use_inertial_stiffness_matrix_K = true; // default true. Can be disabled globally, for testing or optimization
//bool ChLoadAddedMass::use_gyroscopic_torque = true;           // default true. Can be disabled globally, for testing or optimization


ChLoadAddedMass::ChLoadAddedMass(std::shared_ptr<ChBody> body,  ///< object to apply additional inertia to
	const ChVector<>& m_offset,      ///< offset of the center of mass, in body coordinate system
	const double m_mass/*,             ///< added mass [kg]
	const ChVector<>& m_IXX,         ///< added diag. inertia values Ixx, Iyy, Izz (in body coordinate system, centered in body)
	const ChVector<>& m_IXY*/)
	: ChLoadCustom(body), c_m(m_offset), mass(m_mass)
{
	//this->SetInertiaXX(m_IXX);
	//this->SetInertiaXY(m_IXY);
}

// The inertia tensor functions

//void ChLoadAddedMass::SetInertia(const ChMatrix33<>& newXInertia) {
//	I = newXInertia;
//}

//void ChLoadAddedMass::SetInertiaXX(const ChVector<>& iner) {
//	I(0, 0) = iner.x();
//	I(1, 1) = iner.y();
//	I(2, 2) = iner.z();
//}

//void ChLoadAddedMass::SetInertiaXY(const ChVector<>& iner) {
//	I(0, 1) = iner.x();
//	I(0, 2) = iner.y();
//	I(1, 2) = iner.z();
//	I(1, 0) = iner.x();
//	I(2, 0) = iner.y();
//	I(2, 1) = iner.z();
//}

//ChVector<> ChLoadAddedMass::GetInertiaXX() const {
//	ChVector<> iner;
//	iner.x() = I(0, 0);
//	iner.y() = I(1, 1);
//	iner.z() = I(2, 2);
//	return iner;
//}

//ChVector<> ChLoadAddedMass::GetInertiaXY() const {
//	ChVector<> iner;
//	iner.x() = I(0, 1);
//	iner.y() = I(0, 2);
//	iner.z() = I(1, 2);
//	return iner;
//}

void ChLoadAddedMass::ComputeQ(ChState* state_x, ChStateDelta* state_w) { // state_x is position, state_w is velocity?
	//auto mbody = std::dynamic_pointer_cast<ChBody>(this->loadable);
	//if (!mbody->Variables().IsActive())
	//	return;

	//// fetch speeds/pos/accel as 3d vectors for convenience
	//ChVector<> v_x = state_w->segment(0, 3); // abs. 
	//ChVector<> v_w = state_w->segment(3, 3); // local 

	///* // NO ACCELERATION PROPORTIONAL TERM ADDED HERE! Can use LoadIntLoadResidual_Mv if needed.
	//ChVector<> Ma_x = this->mass * (a_x + chrono::Vcross(a_w, this->c_m));
	//ChVector<> Ma_w = this->mass * chrono::Vcross(this->c_m, a_x) + this->I * a_w;
	//*/

	//// Terms of inertial quadratic type (centrifugal, gyroscopic)
	//ChVector<> quadratic_x;
	//ChVector<> quadratic_w;
	//quadratic_x = this->mass * chrono::Vcross(v_w, chrono::Vcross(v_w, this->c_m)); // centrifugal: m*(w X w X c_m)
	///*if (this->use_gyroscopic_torque)*/if(true)
	//	quadratic_w = chrono::Vcross(v_w, this->I * v_w); // gyroscopical: w X J*w
	//else
	//	quadratic_w = VNULL;

	//load_Q.segment(0, 3) = -(quadratic_x).eigen(); // sign: negative, as Q goes in RHS
	//load_Q.segment(3, 3) = -(quadratic_w).eigen(); // sign: negative, as Q goes in RHS
}


void ChLoadAddedMass::ComputeJacobian(ChState* state_x,       ///< state position to evaluate jacobians
	ChStateDelta* state_w,  ///< state speed to evaluate jacobians
	ChMatrixRef mK,         ///< result dQ/dx
	ChMatrixRef mR,         ///< result dQ/dv
	ChMatrixRef mM          ///< result dQ/da
) {
	// fetch speeds as 3d vectors for convenience
	//ChVector<> v_x = state_w->segment(0, 3); // abs. 
	//ChVector<> v_w = state_w->segment(3, 3); // local 
	// (note: accelerations should be fetched from a "state_acc" like we did for speeds with state_w, 
	// but acc. are not available in ComputeQ inputs... maybe in future we should change the ChLoad 
	// class to support also this. For this special case, it works also with the following trick, but 
	// would fail the default numerical differentiation in ComputeJacobian() for the M=dQ/da matrix, so 
	// we override ComputeJacobian and provide the analytical jacobians)
	//auto mbody = std::dynamic_pointer_cast<ChBody>(this->loadable);
	//ChVector<> a_x = mbody->GetA().transpose() * mbody->GetPos_dtdt(); // local 
	//ChVector<> a_w = mbody->GetWacc_loc(); // local 

	//ChStarMatrix33<> wtilde(v_w);  // [w~]
	//ChStarMatrix33<> ctilde(c_m);  // [c~]

	// Analytic expression of inertial load jacobians.
	// Note signs: positive as they go in LHS. 

	// M mass matrix terms (6x6, split in four 3x3 blocks for convenience)
	//jacobians->M.setZero();
	//jacobians->M.block(0, 0, 3, 3).diagonal().setConstant(this->mass);
	//jacobians->M.block(0, 3, 3, 3) = this->mass * chrono::ChStarMatrix33<>(-this->c_m);
	//jacobians->M.block(3, 0, 3, 3) = this->mass * chrono::ChStarMatrix33<>(this->c_m);
	//jacobians->M.block(3, 3, 3, 3) = this->I;

	//set mass matrix here
	jacobians->M(0, 0) = 71.57346 * 1000;
	jacobians->M(1, 1) = 71.57342 * 1000;
	jacobians->M(2, 2) = 130.8768 * 1000;
	jacobians->M(3, 3) = 286.2663 * 1000;
	jacobians->M(4, 4) = 286.2663 * 1000;
	jacobians->M(5, 5) = 6.817555E-8 * 1000;
	jacobians->M(3, 1) = -143.14 * 1000;
	jacobians->M(4, 0) = 143.1401 * 1000;
	jacobians->M(1, 3) = -143.1401 * 1000;
	jacobians->M(0, 4) = 143.1402 * 1000;

	// R gyroscopic damping matrix terms (6x6, split in 3x3 blocks for convenience)
	jacobians->R.setZero();
	///*if (this->use_inertial_damping_matrix_R)*/if(true) {
	//	//  Ri = [0, - m*[w~][c~] - m*[([w~]*c)~]  ; 0 , [w~][I] - [([I]*w)~]  ]
	//	jacobians->R.block(0, 3, 3, 3) = -this->mass * (wtilde * ctilde + ChStarMatrix33<>(wtilde * c_m));
	//	jacobians->R.block(3, 3, 3, 3) = wtilde * I - ChStarMatrix33<>(I * v_w);
	//}

	// K inertial stiffness matrix terms (6x6, split in 3x3 blocks for convenience)
	jacobians->K.setZero();
	///*if (this->use_inertial_stiffness_matrix_K)*/if(true) {
	//	ChStarMatrix33<> atilde(a_w);  // [a~]
	//	// Ki_al = [0, -m*[([a~]c)~] -m*[([w~][w~]c)~] ; 0, m*[c~][xpp~] ]
	//	jacobians->K.block(0, 3, 3, 3) = -this->mass * ChStarMatrix33<>(atilde * c_m) - this->mass * ChStarMatrix33<>(wtilde * (wtilde * c_m));
	//	jacobians->K.block(3, 3, 3, 3) = this->mass * ctilde * ChStarMatrix33<>(a_x);
	//}
}

// The default base implementation in ChLoadCustom could suffice, but here reimplement it in sake of higher speed
// because we can exploiti the sparsity of the formulas.
void ChLoadAddedMass::LoadIntLoadResidual_Mv(ChVectorDynamic<>& R, const ChVectorDynamic<>& w, const double c) {
	if (!this->jacobians)
		return;

	if (!loadable->IsSubBlockActive(0))
		return;

	// fetch w as a contiguous vector
	ChVector<> a_x = w.segment(loadable->GetSubBlockOffset(0), 3);
	ChVector<> a_w = w.segment(loadable->GetSubBlockOffset(0) + 3, 3);

	// R+=c*M*a  
	R.segment(loadable->GetSubBlockOffset(0), 3) += c * (this->mass * (a_x + chrono::Vcross(a_w, this->c_m))).eigen();
	R.segment(loadable->GetSubBlockOffset(0) + 3, 3) += c * (this->mass * chrono::Vcross(this->c_m, a_x) + this->I * a_w).eigen();
}