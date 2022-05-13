#include <cstdio>

#include "chrono/solver/ChSolverPMINRES.h"
#include "chrono/assets/ChPointPointDrawing.h"

#include "chrono/physics/ChSystemNSC.h"
#include "chrono/solver/ChIterativeSolverLS.h"
#include "chrono/timestepper/ChTimestepper.h"
#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChBodyEasy.h"

#include "chrono_irrlicht/ChIrrApp.h"
#include "chrono_irrlicht/ChIrrMeshTools.h"

#include "chrono/fea/ChMeshFileLoader.h"
#include "chrono/physics/ChForce.h"
#include "chrono/physics/ChLoadContainer.h"
#include "chrono/physics/ChLoadsBody.h"

#include "H5Cpp.h"

using namespace chrono;
using namespace chrono::irrlicht;
using namespace chrono::fea;

// =============================================================================
class H5FileInfo {
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
	std::string bodyNum;
	void readH5Data();

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
	//double GetExcitationReal(int m, int n, int w) const;
	//double GetExcitationImag(int m, int n, int w) const;
	double GetRIRFdt() const;
	std::vector<double> GetRIRFTimeVector() const;
	double GetNumFreqs() const;
};

// =============================================================================
class HydroInputs {
public:
	HydroInputs();
	double regular_wave_amplitude;
	double regular_wave_omega;
private:
};

// =============================================================================
class HydroForces;
class ForceTorqueFunc : public ChFunction {
private:
	HydroForces* base;
	int index;

public:
	ForceTorqueFunc(HydroForces* b, int i);
	virtual ForceTorqueFunc* Clone() const override;
	virtual double 	Get_y(double x) const override;
	void SetBase(HydroForces* b);
	void SetIndex(int i);
};
// =============================================================================
class HydroForces {
private:
	std::shared_ptr<ChBody> body;
	H5FileInfo file_info;
	HydroInputs hydro_inputs;
	ChVectorN<double, 6> equilibrium;
	ForceTorqueFunc forces[6];
	std::shared_ptr<ForceTorqueFunc> force_ptrs[6];
	// ChVectorN<double, 6> active_dofs;
	ChVectorN<double, 6> force_hydrostatic;
	ChVectorN<double, 6> force_radiation_damping;
	// ChVectorN<double, 6> forceRadiationDampingFreq;
	ChVectorN<double, 6> force_excitation_freq;
	double wave_amplitude;
	double wave_omega;
	double wave_omega_delta;
	double freq_index_des;
	int freq_index_floor;
	double freq_interp_val;
	ChVectorN<double, 6> excitation_force_mag;
	ChVectorN<double, 6> excitation_force_phase;
	std::vector<ChVectorN<double, 6>> velocity_history;
	double previous_time;
	double previous_time_rirf;
	double previous_time_ex;
	std::vector<double> rirf_time_vector;
	int offset;
	//double current_time;
	std::shared_ptr<ChForce> chrono_force;
	std::shared_ptr<ChForce> chrono_torque;
public:
	HydroForces();
	HydroForces(H5FileInfo& h5_file_info, std::shared_ptr<ChBody> object, HydroInputs users_hydro_inputs);

	// ensures these member functions are not given a default definition by compiler
	// these features don't make sense to have, and break things when they exist
	HydroForces(const HydroForces& other) = delete;
	HydroForces operator = (const HydroForces& rhs) = delete;

	ChVectorN<double, 6> ComputeForceHydrostatics();
	ChVectorN<double, 6> ComputeForceRadiationDampingConv();
	// ChVectorN<double, 6> ComputeForceRadiationDampingFreq();
	ChVectorN<double, 6> ComputeForceExcitationRegularFreq();
	double coordinateFunc(int i);
	void SetForce();
	void SetTorque();
};

// =============================================================================
class ChLoadAddedMass : public ChLoadCustom {
public:
	ChLoadAddedMass(std::shared_ptr<ChBody> body,  ///< object to apply additional inertia to
		const H5FileInfo& file );                ///< h5 file to initialize added mass with

	/// "Virtual" copy constructor (covariant return type).
	virtual ChLoadAddedMass* Clone() const override { return new ChLoadAddedMass(*this); }

	/// Compute Q, the generalized load.
	/// In this case, it computes the quadratic (centrifugal, gyroscopic) terms.
	/// Signs are negative as Q assumed at right hand side, so Q= -Fgyro -Fcentrifugal
	/// Called automatically at each Update().
	/// The M*a term is not added: to this end one could use LoadIntLoadResidual_Mv afterward.
	virtual void ComputeQ(ChState* state_x,      ///< state position to evaluate Q
		ChStateDelta* state_w  ///< state speed to evaluate Q
	) override {}

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
	ChMatrixDynamic<double> inf_added_mass_J;       ///< added mass at infinite frequency in global coordinates

	virtual bool IsStiff() override { return true; } // this to force the use of the inertial M, R and K matrices

};
// =============================================================================
class LoadAllHydroForces {
private:
	H5FileInfo sys_file_info;     /// < object to read h5 file info
	HydroForces hydro_force;                     /// < object for linear restoring force
	HydroInputs users_hydro_inputs;
	//std::shared_ptr<BuoyancyForce> buoyancy_force;
	std::shared_ptr<ChLoadContainer> my_loadcontainer;
	std::shared_ptr<ChLoadAddedMass> my_loadbodyinertia;
public:
	LoadAllHydroForces(std::shared_ptr<ChBody> object, std::string file, std::string body_name, HydroInputs users_hydro_inputs);
	//~LoadAllHydroForces();
};
