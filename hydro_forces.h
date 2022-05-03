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
	std::vector<double> timesteps;
	hsize_t freq_dims[3];
	std::vector<double> freq_list;
	double omega_min;
	double omega_max;
	double rho;
	double g;
	double disp_vol;
	std::string h5_file_name;
	std::string bodyNum;
	void read_data();

public:
	H5FileInfo();
	H5FileInfo(std::string file, std::string bodyName); //TODO: systemName
	~H5FileInfo();
	ChMatrixDynamic<double> get_lin_matrix() const;
	ChMatrixDynamic<double> get_inf_added_mass_matrix() const;
	ChVector<> get_equil_cg() const;
	ChVector<> get_equil_cb() const;
	double get_rho() const;
	double get_g() const;
	double get_disp_vol() const;
	double get_rirf_val(int i, int n, int m) const;
	double get_excitation_mag(int m, int n, int w) const;
	double get_excitation_phase(int m, int n, int w) const;
	double get_omega_min() const;
	double get_omega_max() const;
	double get_domega() const;
	//double get_excitation_re(int m, int n, int w) const;
	//double get_excitation_im(int m, int n, int w) const;
	int get_rirf_dims(int i) const;
	double get_delta_t() const;
	std::vector<double> get_times() const;
	double get_num_freqs() const;
};

// =============================================================================
class HydroInputs {
public:
	HydroInputs();
	//~HydroInputs();
	double regularWaveAmplitude;
	double regularWaveOmega;
	int freqIndex;
	//double regularWavePeriod;
	//double regularWaveOmega;
	//double get_regular_wave_omega(double regularWavePeriod);
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
	H5FileInfo fileInfo;
	HydroInputs hydroInputs;
	ChVectorN<double, 6> equil;
	ForceTorqueFunc forces[6];
	std::shared_ptr<ForceTorqueFunc> force_ptrs[6];
	ChVectorN<double, 6> activeDofs;
	ChVectorN<double, 6> currentForce;
	ChVectorN<double, 6> forceRadiationDamping;
	ChVectorN<double, 6> forceRadiationDampingFreq;
	ChVectorN<double, 6> forceExcitation;
	double waveAmplitude;
	double waveOmega;
	double domega;
	double freqIndexDes;
	int freqIndexFloor;
	double freqInterpVal;
	double forceExcitationMag;
	double forceExcitationPhase;
	std::vector<ChVectorN<double, 6>> velocityHistory;
	std::vector<double> timeSteps;
	int offset;
	double currentTime;
	double prevTime;
	double prevTimeIRF;
	double prevTimeEx;
	std::shared_ptr<ChForce> chronoForce;
	std::shared_ptr<ChForce> chronoTorque;
public:
	HydroForces();
	HydroForces(H5FileInfo& sysH5FileInfo, std::shared_ptr<ChBody> object, HydroInputs userHydroInputs);

	// ensures these member functions are not given a default definition by compiler
	// these features don't make sense to have, and break things when they exist
	HydroForces(const HydroForces& other) = delete;
	HydroForces operator = (const HydroForces& rhs) = delete;

	ChVectorN<double, 6> fHydrostaticStiffness();
	ChVectorN<double, 6> fRadDamping();
	ChVectorN<double, 6> fRadDampingFreq();
	ChVectorN<double, 6> fExcitationRegularFreq();
	double coordinateFunc(int i);
	void SetForce();
	void SetTorque();
};
// TDO: bring this back!
//// =============================================================================
//class BuoyancyForce {
//private:
//	H5FileInfo fileInfo;
//	double buoyancy;
//	ChFunction_Const fc;
//	std::shared_ptr<ChFunction_Const> fc_ptr;
//	ChForce force;
//	std::shared_ptr<ChForce> force_ptr;
//
//public:
//	BuoyancyForce(H5FileInfo& info);
//	std::shared_ptr<ChForce> getForce_ptr();
//};
// =============================================================================

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
	H5FileInfo sysFileInfo;     /// < object to read h5 file info
	HydroForces hydro_force;                     /// < object for linear restoring force
	HydroInputs userHydroInputs;
	//std::shared_ptr<BuoyancyForce> buoyancy_force;
	std::shared_ptr<ChLoadContainer> my_loadcontainer;
	std::shared_ptr<ChLoadAddedMass> my_loadbodyinertia;
public:
	LoadAllHydroForces(std::shared_ptr<ChBody> object, std::string file, HydroInputs userHydroInputs);
};
