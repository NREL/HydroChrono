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


