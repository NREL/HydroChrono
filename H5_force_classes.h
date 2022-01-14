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

#include "chrono/fea/ChElementAddedMass.h"
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
	ChVector<double> cg;
	ChVector<double> cb;
	double rho;
	double g;
	double disp_vol;
	std::string h5_file_name;
	std::string bodyNum;
	void read_data(); 

public:
	BodyFileInfo();
	BodyFileInfo(std::string file, std::string bodyName);
	ChMatrixDynamic<double> get_lin_matrix() const;
	ChVector<> get_equil_cg() const;
	ChVector<> get_equil_cb() const;
	double get_rho() const;
	double get_g() const;
	double get_disp_vol() const;
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

public:
	LinRestorForce();
	LinRestorForce(BodyFileInfo& lin, std::shared_ptr<ChBody> object);

	// ensures these member functions are not given a default definition by compiler
	// these features don't make sense to have, and break things when they exist
	LinRestorForce(const LinRestorForce& other) = delete;
	LinRestorForce operator = (const LinRestorForce& rhs) = delete;

	ChVectorN<double, 6> Get_p() const;

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

