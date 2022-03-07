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

//#include "chrono/fea/ChElementAddedMass.h"
#include "chrono/fea/ChNodeFEAxyz.h"
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
	double* K_matrix;
	hsize_t K_dims[3];
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
	ChVector<> get_equil_cg() const;
	ChVector<> get_equil_cb() const;
	double get_rho() const;
	double get_g() const;
	double get_disp_vol() const;
	double get_impulse_resp(int i, int n, int m) const;
	int get_K_dims(int i) const;
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

class LinRestorForce {
private:
	std::shared_ptr<ChBody> bobber;
	BodyFileInfo fileInfo;
	ChVectorN<double, 6> equil;
	ForceTorqueFunc forces[6];
	std::shared_ptr<ForceTorqueFunc> force_ptrs[6];
	ChVectorN<double, 6> currentForce;
	double prevTime;

public:
	LinRestorForce();
	LinRestorForce(BodyFileInfo& lin, std::shared_ptr<ChBody> object);

	// ensures these member functions are not given a default definition by compiler
	// these features don't make sense to have, and break things when they exist
	LinRestorForce(const LinRestorForce& other) = delete;
	LinRestorForce operator = (const LinRestorForce& rhs) = delete;

	ChVectorN<double, 6> matrixMult();

	double coordinateFunc(int i);
	void SetForce(std::shared_ptr<ChForce> force);
	void SetTorque(std::shared_ptr<ChForce> torque);
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

public:
	ImpulseResponseForce();
	ImpulseResponseForce(BodyFileInfo& file, std::shared_ptr<ChBody> object);
	// ensures these member functions are not given a default definition by compiler
	// these features don't make sense to have, and break things when they exist
	ImpulseResponseForce(const ImpulseResponseForce& other) = delete;
	ImpulseResponseForce operator = (const ImpulseResponseForce& rhs) = delete;

	ChVectorN<double, 6> convolutionIntegral(); 

	double coordinateFunc(int i);
	void SetForce(std::shared_ptr<ChForce> force);
	void SetTorque(std::shared_ptr<ChForce> torque);
};

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

