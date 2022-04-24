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
class BodyFileInfo {
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
	double rho;
	double g;
	double disp_vol;
	std::string h5_file_name;
	std::string bodyNum;
	void read_data();

public:
	BodyFileInfo();
	BodyFileInfo(std::string file, std::string bodyName);
	~BodyFileInfo();
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
	double get_excitation_re(int m, int n, int w) const;
	double get_excitation_im(int m, int n, int w) const;
	int get_rirf_dims(int i) const;
	double get_delta_t() const;
	std::vector<double> get_times() const;
};
// =============================================================================
class LinRestorForce;
class ForceTorqueFunc : public ChFunction {
private:
	LinRestorForce* base;
	int index;

public:
	ForceTorqueFunc(LinRestorForce* b, int i);
	virtual ForceTorqueFunc* Clone() const override;
	virtual double 	Get_y(double x) const override;
	void SetBase(LinRestorForce* b);
	void SetIndex(int i);
};
// =============================================================================
class LinRestorForce {
private:
	std::shared_ptr<ChBody> bobber;
	BodyFileInfo fileInfo;
	ChVectorN<double, 6> equil;
	ForceTorqueFunc forces[6];
	std::shared_ptr<ForceTorqueFunc> force_ptrs[6];
	ChVectorN<double, 6> currentForce;
	double prevTime;
	std::shared_ptr<ChForce> placeholder;
	std::shared_ptr<ChForce> placeholdertorque;
public:
	LinRestorForce();
	LinRestorForce(BodyFileInfo& lin, std::shared_ptr<ChBody> object);

	// ensures these member functions are not given a default definition by compiler
	// these features don't make sense to have, and break things when they exist
	LinRestorForce(const LinRestorForce& other) = delete;
	LinRestorForce operator = (const LinRestorForce& rhs) = delete;

	ChVectorN<double, 6> matrixMult();

	double coordinateFunc(int i);
	void SetForce();
	void SetTorque();
};
// =============================================================================
class BuoyancyForce {
private:
	BodyFileInfo fileInfo;
	double bf;
	ChFunction_Const fc;
	std::shared_ptr<ChFunction_Const> fc_ptr;
	ChForce force;
	std::shared_ptr<ChForce> force_ptr;

public:
	BuoyancyForce(BodyFileInfo& info);
	std::shared_ptr<ChForce> getForce_ptr();
};
// =============================================================================
class ImpulseResponseForce;
class IRF_func : public ChFunction {
private:
	ImpulseResponseForce* base;
	int index;

public:
	IRF_func(ImpulseResponseForce* b, int i);
	virtual IRF_func* Clone() const override;
	virtual double 	Get_y(double x) const override;
	void SetBase(ImpulseResponseForce* b);
	void SetIndex(int i);
};
// =============================================================================
class ImpulseResponseForce {
private:
	BodyFileInfo fileInfo;
	std::shared_ptr<ChBody> body;
	std::vector<ChVectorN<double, 6>> velHistory;
	ChVectorN<double, 6> currentForce;
	std::vector<double> timeSteps;
	IRF_func forces[6];
	std::shared_ptr<IRF_func> force_ptrs[6];
	int offset;
	double prevTime;
	std::shared_ptr<ChForce> placeholder;
	std::shared_ptr<ChForce> placeholdertorque;
public:
	ImpulseResponseForce();
	ImpulseResponseForce(BodyFileInfo& file, std::shared_ptr<ChBody> object);
	// ensures these member functions are not given a default definition by compiler
	// these features don't make sense to have, and break things when they exist
	ImpulseResponseForce(const ImpulseResponseForce& other) = delete;
	ImpulseResponseForce operator = (const ImpulseResponseForce& rhs) = delete;

	ChVectorN<double, 6> convolutionIntegral();

	double coordinateFunc(int i);
	void SetForce();
	void SetTorque();
};
// =============================================================================
class ChLoadAddedMass : public ChLoadCustom {
public:
	ChLoadAddedMass(std::shared_ptr<ChBody> body,  ///< object to apply additional inertia to
		const BodyFileInfo& file );                ///< h5 file to initialize added mass with

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
	BodyFileInfo file_info;     /// < object to read h5 file info
	LinRestorForce lin_restor_force_2;                     /// < object for linear restoring force
	ImpulseResponseForce irf;                              /// < object for impulse restoring force

	// declare some forces to be initialized in lin_restor_force_2 to be applied to to body later
	//std::shared_ptr<ChForce> force;
	//std::shared_ptr<ChForce> torque;
	//std::shared_ptr<ChForce> force2;
	//std::shared_ptr<ChForce> torque2;
	std::shared_ptr<BuoyancyForce> fb;
	std::shared_ptr<ChLoadContainer> my_loadcontainer;
	std::shared_ptr<ChLoadAddedMass> my_loadbodyinertia;
public:
	LoadAllHydroForces(std::shared_ptr<ChBody> object, std::string file);
};
