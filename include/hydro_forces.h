#ifndef _HYDRO_FORCES_
#define _HYDRO_FORCES_

#ifndef SWIG
#include <filesystem> // includes file checking features for h5 class
#include "H5Cpp.h"
#include "chrono/physics/ChSystemNSC.h" // includes ChSystem which includes many std libs and chrono core files
#include "chrono/physics/ChLoad.h"
#include "chrono/physics/ChLoadContainer.h" // includes ChLoad and ChPhysicsItem
#include "chrono/physics/ChLoadsBody.h" // includes ChFunction, ChBody, and ChLoad
#endif

#include "swig_test.h"

#ifdef hydro_lib_EXPORTS
#define hydro_lib_API __declspec(dllexport)
#else
#define hydro_lib_API 
#endif

using chrono::ChMatrixDynamic;
using chrono::ChBody;
using chrono::ChForce;
using chrono::ChLoadCustomMultiple;
using chrono::ChLoadable;

// =============================================================================
class hydro_lib_API H5FileInfo {
public:
	bool printed = false;
	H5FileInfo();
	H5FileInfo(std::string file, std::string body_name);
	H5FileInfo(const H5FileInfo& old);
	H5FileInfo& operator = (const H5FileInfo& rhs);
	void InitScalar(H5::H5File& file, std::string data_name, double& var);
	void Init1D(H5::H5File& file, std::string data_name, std::vector<double>& var);
	void Init2D(H5::H5File& file, std::string data_name, ChMatrixDynamic<double>& var); 
	void Init3D(H5::H5File& file, std::string data_name, std::vector<double>& var, std::vector<int>& dims);
	~H5FileInfo();

	ChMatrixDynamic<double> GetInfAddedMassMatrix() const;
	double GetHydrostaticStiffness(int i, int j) const;
	double GetRIRFval(int i, int n, int m) const;
	int GetRIRFDims(int i) const;
	std::vector<double> GetRIRFTimeVector() const; // TODO
	double GetExcitationMagValue(int m, int n, int w) const;
	double GetExcitationMagInterp(int i, int j, double freq_index_des) const;
	double GetOmegaDelta() const;
	double GetOmegaMax() const;
	double GetExcitationPhaseValue(int m, int n, int w) const;
	double GetExcitationPhaseInterp(int i, int j, double freq_index_des) const;
	double GetNumFreqs() const;

	std::vector<double> cg;
	std::vector<double> cb;
	const double& rho = _rho;
	const double& g = _g;
	const double& disp_vol = _disp_vol;
	//const double& rirf_timestep = _rirf_timestep;
	int bodyNum;
private:
	double _rho;
	double _g;
	double _disp_vol;
	double _rirf_timestep;
	std::vector<double> freq_list;
	ChMatrixDynamic<double> lin_matrix;
	ChMatrixDynamic<double> inf_added_mass;
	std::vector<double> rirf_matrix;
	std::vector<int> rirf_dims;
	std::vector<double> rirf_time_vector;
	std::vector<double> radiation_damping_matrix; // TODO check about names
	std::vector<int> Bw_dims; // TODO check with dave on name for dimensions of radiation damping matrix
	std::vector<double> excitation_mag_matrix;
	std::vector<int> excitation_mag_dims;
	std::vector<double> excitation_re_matrix;
	std::vector<int> re_dims;
	std::vector<double> excitation_im_matrix;
	std::vector<int> im_dims;
	std::vector<double> excitation_phase_matrix;
	std::vector<int> excitation_phase_dims;
	std::string h5_file_name;
	std::string bodyName;
	void readH5Data();
};


// =============================================================================
class hydro_lib_API ForceFunc6d;
class HydroInputs;
class hydro_lib_API TestHydro;

class hydro_lib_API ComponentFunc : public chrono::ChFunction {
public:
	ComponentFunc();
	ComponentFunc(const ComponentFunc& old);
	ComponentFunc(ForceFunc6d* b, int i);
	virtual ComponentFunc* Clone() const override;
	virtual double 	Get_y(double x) const override;
private:
	ForceFunc6d* base; 
	int index;
};

// =============================================================================
// ForceFunc6d organizes the functional (time dependent) forces in each DoF (6 total) for a body
class hydro_lib_API ForceFunc6d {
public:
	ForceFunc6d();
	ForceFunc6d(std::shared_ptr<ChBody> object, TestHydro* all_hydro_forces_user);
	ForceFunc6d(const ForceFunc6d& old);
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

class hydro_lib_API ChLoadAddedMass;


// =============================================================================
class hydro_lib_API ChLoadAddedMass : public ChLoadCustomMultiple {
public:
	ChLoadAddedMass(const std::vector<H5FileInfo>& file,   ///< h5 file to initialize added mass with
					std::vector<std::shared_ptr<ChBody>>& bodies  ///< objects to apply additional inertia to
					);     
//	/// "Virtual" copy constructor (covariant return type).
	virtual ChLoadAddedMass* Clone() const override { return new ChLoadAddedMass(*this); }
//
//	/// Compute Q, the generalized load.
//	/// In this case, it computes the quadratic (centrifugal, gyroscopic) terms.
//	/// Signs are negative as Q assumed at right hand side, so Q= -Fgyro -Fcentrifugal
//	/// Called automatically at each Update().
//	/// The M*a term is not added: to this end one could use LoadIntLoadResidual_Mv afterward.
	virtual void ComputeQ(chrono::ChState* state_x,      ///< state position to evaluate Q
		chrono::ChStateDelta* state_w  ///< state speed to evaluate Q
	) override {}

	/// For efficiency reasons, do not let the parent class do automatic differentiation
	/// to compute the R, K matrices. Use analytic expressions instead. For example, R is
	/// the well known gyroscopic damping matrix. Also, compute the M matrix.
	virtual void ComputeJacobian(chrono::ChState* state_x,       ///< state position to evaluate jacobians
		chrono::ChStateDelta* state_w,  ///< state speed to evaluate jacobians
		chrono::ChMatrixRef mK,         ///< result -dQ/dx
		chrono::ChMatrixRef mR,         ///< result -dQ/dv
		chrono::ChMatrixRef mM          ///< result -dQ/da
	) override;

	/// Just for efficiency, override the default LoadIntLoadResidual_Mv, because we can do this in a simplified way.
	virtual void LoadIntLoadResidual_Mv(chrono::ChVectorDynamic<>& R,           ///< result: the R residual, R += c*M*w
		const chrono::ChVectorDynamic<>& w,     ///< the w vector
		const double c) override;       ///< a scaling factor
	void AssembleSystemAddedMassMat();
private:
	ChMatrixDynamic<double> infinite_added_mass;       ///< added mass at infinite frequency in global coordinates
	int nBodies;
	std::vector<H5FileInfo> h5_body_data;
	virtual bool IsStiff() override { return true; } // this to force the use of the inertial M, R and K matrices
};
#endif