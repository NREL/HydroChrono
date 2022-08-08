#include <cstdio>

#include "chrono/solver/ChSolverPMINRES.h"
#include "chrono/solver/ChIterativeSolverLS.h"
#include "chrono/timestepper/ChTimestepper.h"

#include "chrono/physics/ChForce.h"
#include "chrono/physics/ChLoadContainer.h"
#include "chrono/physics/ChLoadsBody.h"
#include "chrono/physics/ChLoad.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChBodyEasy.h"

#include "chrono/fea/ChMeshFileLoader.h"

#include "chrono/assets/ChPointPointDrawing.h"
#include "chrono_irrlicht/ChIrrApp.h"
#include "chrono_irrlicht/ChIrrMeshTools.h"

#include "H5Cpp.h"

#include <filesystem>

using namespace chrono;
using namespace chrono::irrlicht;
using namespace chrono::fea;

// =============================================================================
class H5FileInfo { // TODO cut simple setter/getter functions to trim this huge class down
public:
	H5FileInfo();
	H5FileInfo(std::string file, std::string body_name);
	~H5FileInfo();
	ChMatrixDynamic<double> GetHydrostaticStiffnessMatrix() const;
	ChMatrixDynamic<double> GetInfAddedMassMatrix() const;
	ChVector<> GetEquilibriumCoG() const;
	ChVector<> GetEquilibriumCoB() const;
	double GetRho() const;
	double GetGravity() const;
	double GetDisplacementVolume() const;
	double GetRIRFval(int i, int n, int m) const;
	int GetRIRFDims(int i) const;
	double GetExcitationMagValue(int m, int n, int w) const;
	double GetExcitationMagInterp(int i, int j, double freq_index_des) const;
	double GetExcitationPhaseValue(int m, int n, int w) const;
	double GetExcitationPhaseInterp(int i, int j, double freq_index_des) const;
	double GetOmegaMin() const;
	double GetOmegaMax() const;
	double GetOmegaDelta() const;
	double GetRIRFdt() const;
	std::vector<double> GetRIRFTimeVector() const;
	double GetNumFreqs() const;
	int bodyNum;
private:
	ChMatrixDynamic<double> lin_matrix;
	ChMatrixDynamic<double> inf_added_mass;
	double* rirf_matrix;
	hsize_t rirf_dims[3];
	double* radiation_damping_matrix;
	hsize_t radiation_damping_dims[3];
	double* excitation_mag_matrix;
	hsize_t excitation_mag_dims[3];
	double* excitation_phase_matrix;
	hsize_t excitation_phase_dims[3];
	double* excitation_re_matrix;
	hsize_t excitation_re_dims[3];
	double* excitation_im_matrix;
	hsize_t excitation_im_dims[3];
	ChVector<double> cg;
	ChVector<double> cb;
	std::vector<double> rirf_time_vector;
	hsize_t freq_dims[3];
	std::vector<double> freq_list;
	double omega_min;
	double omega_max;
	double rho;
	double g;
	double disp_vol;
	std::string h5_file_name;
	std::string bodyName;
	void readH5Data();
};

// =============================================================================
class HydroInputs {
public:
	HydroInputs();
	double SetRegularWaveAmplitude(double val) {
		regular_wave_amplitude = val;
		return regular_wave_amplitude;
	}
	double GetRegularWaveAmplitude() const { return regular_wave_amplitude; }
	double SetRegularWaveOmega(double val) {
		regular_wave_omega = val;
		return regular_wave_omega;
	}
	double GetRegularWaveOmega() const { return regular_wave_omega; }
	
private:
	double regular_wave_amplitude;
	double regular_wave_omega;
};

// =============================================================================
class ForceFunc6d;

class TestHydro;

class ComponentFunc : public ChFunction {
public:
	ComponentFunc(ForceFunc6d* b, int i);
	virtual ComponentFunc* Clone() const override;
	virtual double 	Get_y(double x) const override;
private:
	ForceFunc6d* base; 
	int index;
};

// =============================================================================
// ForceFunc6d organizes the functional (time dependent) forces in each DoF (6 total) for a body
class ForceFunc6d {
public:
	ForceFunc6d();
	ForceFunc6d(std::shared_ptr<ChBody> object, TestHydro* all_hydro_forces_user);
	//ForceFunc6d(const ForceFunc6d& other) = delete;
	ForceFunc6d operator = (const ForceFunc6d& rhs) = delete;
	double coordinateFunc(int i);
private:
	void SetForce();
	void SetTorque();
	std::shared_ptr<ChBody> body;
	int b_num;
	ComponentFunc forces[6];
	std::shared_ptr<ComponentFunc> force_ptrs[6];
	std::shared_ptr<ChForce> chrono_force;
	std::shared_ptr<ChForce> chrono_torque;
	TestHydro* all_hydro_forces;
};

class TestHydro {
public:
	TestHydro();
	TestHydro(std::vector<std::shared_ptr<ChBody>> user_bodies, std::string h5_file_name, HydroInputs users_hydro_inputs);
	TestHydro(const TestHydro& other) = delete;
	TestHydro operator = (const TestHydro& rhs) = delete;
	std::vector<double> ComputeForceHydrostatics();
	std::vector<double> ComputeForceRadiationDampingConv();
	std::vector<double> ComputeForceExcitation();
	double GetRIRFval(int row, int col, int st);
	double coordinateFunc(int b, int i);
	//ChVectorN<double, 6> ComputeForceExcitationRegularFreq();
private:
	std::vector<std::shared_ptr<ChBody>> bodies;
	std::vector<H5FileInfo> file_info;
	std::vector<ForceFunc6d> force_per_body;
	//HydroInputs hydro_inputs;
	std::vector<double> force_hydrostatic;
	std::vector<double> force_radiation_damping;
	//std::vector<double> force_excitation;
	std::vector<double> total_force;
	int num_bodies;
	std::vector<double> equilibrium;
	double getVelHistoryAllBodies(int step, int c) const;
	double setVelHistory(double val, int step, int b_num, int index);
	//std::vector<double> force_excitation_freq;
	//double wave_amplitude;
	//double wave_omega;
	//double wave_omega_delta;
	//double freq_index_des;
	//int freq_index_floor;
	//double freq_interp_val;
	//ChVectorN<double, 6> excitation_force_mag;
	//ChVectorN<double, 6> excitation_force_phase;
	std::vector<double> velocity_history; // use helper function to access vel_history elments correctly
	double prev_time;
	//double prev_time_rirf;
	//double prev_time_ex;
	std::vector<double> rirf_time_vector; // (should be the same for each body?)
	int offset_rirf;
};

// =============================================================================
//class ChLoadAddedMass : public ChLoadCustomMultiple {
//public:
//	ChLoadAddedMass(const H5FileInfo& file,   ///< h5 file to initialize added mass with
//		std::vector<std::shared_ptr<ChBody>>& bodies  ///< objects to apply additional inertia to
//		);     
//	ChLoadAddedMass(const ChMatrixDynamic<>& addedMassMatrix,   ///< h5 file to initialize added mass with
//		std::vector<std::shared_ptr<ChBody>>& bodies  ///< objects to apply additional inertia to
//	);
//
//	/// "Virtual" copy constructor (covariant return type).
//	virtual ChLoadAddedMass* Clone() const override { return new ChLoadAddedMass(*this); }
//
//	/// Compute Q, the generalized load.
//	/// In this case, it computes the quadratic (centrifugal, gyroscopic) terms.
//	/// Signs are negative as Q assumed at right hand side, so Q= -Fgyro -Fcentrifugal
//	/// Called automatically at each Update().
//	/// The M*a term is not added: to this end one could use LoadIntLoadResidual_Mv afterward.
//	virtual void ComputeQ(ChState* state_x,      ///< state position to evaluate Q
//		ChStateDelta* state_w  ///< state speed to evaluate Q
//	) override {}
//
//	/// For efficiency reasons, do not let the parent class do automatic differentiation
//	/// to compute the R, K matrices. Use analytic expressions instead. For example, R is
//	/// the well known gyroscopic damping matrix. Also, compute the M matrix.
//	virtual void ComputeJacobian(ChState* state_x,       ///< state position to evaluate jacobians
//		ChStateDelta* state_w,  ///< state speed to evaluate jacobians
//		ChMatrixRef mK,         ///< result -dQ/dx
//		ChMatrixRef mR,         ///< result -dQ/dv
//		ChMatrixRef mM          ///< result -dQ/da
//	) override;
//
//	/// Just for efficiency, override the default LoadIntLoadResidual_Mv, because we can do this in a simplified way.
//	virtual void LoadIntLoadResidual_Mv(ChVectorDynamic<>& R,           ///< result: the R residual, R += c*M*w
//		const ChVectorDynamic<>& w,     ///< the w vector
//		const double c) override;       ///< a scaling factor
//private:
//	ChMatrixDynamic<double> inf_added_mass_J;       ///< added mass at infinite frequency in global coordinates
//
//	virtual bool IsStiff() override { return true; } // this to force the use of the inertial M, R and K matrices
//
//};

