// module name is just a useful name for organizing files/namespaces/classes
// each module in a collection is created with separate swig invocations
// when modules need to share information (class derived from another in separate module) it's a bit complex 
// see import module
%module(directors="1") hydro_py_module 

// Python requires the name of module that the base class exists in 
// so that the proxy classes can fully inherit the base class's methods
// desirable to import the header file rather than the interface file and overcome the above warning. 
// probably need to import the module with chbody in it or something
//%import(module="base_module") "base.h"


// use shared libraries is recommended by swig, Due to the complexity of working with shared libraries and multiple modules, it might be a good idea to consult an outside reference. John Levine's "Linkers and Loaders" is highly recommended. 

// Everything in the %{ ... %} block is simply copied verbatim to the resulting wrapper file created by SWIG
// used to include header files and other declarations that are required to make the generated wrapper code compile
// none of this is parsed/interpreted by swig, simply copy pasted to wtapper code
// For enumerations, it is critical that the original enum definition be included somewhere in the interface file 
// (either in a header file or in the %{,%} block). SWIG only translates the enumeration into code needed to add 
// the constants to a scripting language. It needs the original enumeration declaration in order to get the correct 
// enum values as assigned by the C compiler. 
%{
#include "swig_test.h"
#include <vector>
#include "hydro_forces.h" 
%}

%{
#include <typeindex>
#include <cstddef>

#include "chrono/ChConfig.h"

#include "chrono/core/ChApiCE.h"
#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChLink.h"
#include "chrono/physics/ChLinkMotionImposed.h"
#include "chrono/physics/ChLoad.h"
#include "chrono/physics/ChLoadsBody.h"
#include "chrono/physics/ChNodeBase.h"
#include "chrono/physics/ChNodeXYZ.h"
#include "chrono/physics/ChLoadsXYZnode.h"
#include "chrono/physics/ChIndexedNodes.h"

#include "chrono/assets/ChLineShape.h"
#include "chrono/assets/ChPathShape.h"
#include "chrono/assets/ChPointPointShape.h"
#include "chrono/assets/ChSurfaceShape.h"
#include "chrono/assets/ChTriangleMeshShape.h"
#include "chrono/assets/ChEllipsoidShape.h"
#include "chrono/assets/ChVisualMaterial.h"
#include "chrono/assets/ChGlyphs.h"
#include "chrono/assets/ChVisualSystem.h"

#include "chrono/collision/ChCollisionUtils.h"
#include "chrono/collision/ChCollisionSystem.h"
#include "chrono/collision/ChCollisionSystemBullet.h"

#include "chrono/geometry/ChTriangleMesh.h"
#include "chrono/geometry/ChTriangleMeshConnected.h"
#include "chrono/geometry/ChTriangleMeshSoup.h"
#include "chrono/core/ChBezierCurve.h"
#include "Eigen/src/Core/util/Memory.h"
#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono/utils/ChFilters.h"
#include "chrono/utils/ChUtilsCreators.h"
%}

// However, #include statements are ignored unless the -includeall command line option has been supplied.  (out of %{...%} block)
// Now list ANSI C/C++ declarations
// this part is like a preprcessor....

// %inlcude over # includes each file once (dont need include guards)
// To include another file into a SWIG interface, use the %include directive like this:
// this makes more wrapper for pointer.i in this swig 
//%include "pointer.i"
// %feature("valuewrapper") std::vector<(chrono::ChBodyEasyMesh)>;
// %rename("$ignore",$isconstructor,fullname=1) "std::vector<(chrono::ChBodyEasyMesh)>";
%include "std_vector.i"
%include "pyabc.i"
%include "typemaps.i"
%include "std_shared_ptr.i"
%include "exception.i"
// %import(module="pychrono") "chrono_swig/interface/core/ChModuleCore.i"

%shared_ptr(chrono::ChBody)

%shared_ptr(chrono::ChBodyEasySphere)
%shared_ptr(chrono::ChBodyEasyBox)
%shared_ptr(chrono::ChBodyEasyEllipsoid)
%shared_ptr(chrono::ChBodyEasyCylinder)
%shared_ptr(chrono::ChBodyEasyConvexHull)
%shared_ptr(chrono::ChBodyEasyConvexHullAuxRef)
%shared_ptr(chrono::ChBodyEasyMesh)
%shared_ptr(chrono::ChBodyEasyClusterOfSpheres)

%shared_ptr(chrono::ChBodyAuxRef)
%shared_ptr(chrono::cascade::ChCascadeBodyEasy)
%shared_ptr(chrono::cascade::ChCascadeBodyEasyProfile)
%shared_ptr(chrono::ChBodyEasyConvexHullAuxRef)
%shared_ptr(chrono::ChBodyEasyMesh)

// Instantiate templates used by example
%template(vector_double) std::vector<double>;
%template(shared_ChBody)  std::shared_ptr<chrono::ChBody>;
%template(shared_ChBodyEasyMesh) std::shared_ptr<chrono::ChBodyEasyMesh>;
%template(vector_ChBody) std::vector< chrono::ChBody >;
%template(vector_ChBody_ptr) std::vector< std::shared_ptr<chrono::ChBody> >;
%template(vector_ChBodyEasyMesh_ptr) std::vector< std::shared_ptr<chrono::ChBodyEasyMesh> >;
%template(vector_H5FileInfo) std::vector< H5FileInfo >;
%template(vector_ForceFunc6d) std::vector< ForceFunc6d >;

// %inline %{ // add one of these for each body type you wish to upcast to generic ChBody
//   std::vector<std::shared_ptr<chrono::ChBody>> upcast(std::vector<std::shared_ptr<chrono::ChBodyEasyMesh >>& d) {
// 	std::vector<std::shared_ptr<chrono::ChBody>> re(d.size());
// 	std::transform(d.begin(), d.end(), re.begin(), [](const std::shared_ptr<chrono::ChBodyEasyMesh>& p) { return std::static_pointer_cast<chrono::ChBody>(p); });
//     return re;
//   }
// %}
// %inline %{ // add one of these for each body type you wish to upcast to generic ChBody
//   std::shared_ptr<chrono::ChBody> upcast(std::shared_ptr<chrono::ChBodyEasyMesh >& d) {
// 	return d;
//   }
// %}
%inline %{
  std::shared_ptr<chrono::ChBody> upcast(const std::shared_ptr<chrono::ChBodyEasyMesh> d) {
	std::shared_ptr<chrono::ChBody> re = std::static_pointer_cast<chrono::ChBody>(d);
    return re;
  }
%}
// %inline %{ 
//   std::vector<std::shared_ptr<chrono::ChBody *>> upcast(std::vector<std::shared_ptr<chrono::ChBodyEasyMesh *>> d) {
// 	std::vector<std::shared_ptr<chrono::ChBody>> re(d.size());
// 	std::transform(d.begin(), d.end(), re.begin(), [](const std::shared_ptr<chrono::ChBodyEasyMesh>& p) { return std::static_pointer_cast<chrono::ChBody>(p); });
//     return re;
//   }
// %}
/* helper function that helps hydro forces...
std::vector<std::shared_ptr<ChLoadable>> constructorHelper(std::vector<std::shared_ptr<ChBody>>& bodies) {
	std::vector<std::shared_ptr<ChLoadable>> re(bodies.size());
	std::transform(bodies.begin(), bodies.end(), re.begin(), [](const std::shared_ptr<ChBody>& p) { return std::static_pointer_cast<ChLoadable>(p); });
	return re;
}
*/

%feature("director") ChBody;
// %nodefaultctor chrono::ChBodyEasyMesh;
// %import "foo.i" on the other hand would import info from foo.i but not make more wrapper code
// Such information generally includes type declarations (e.g., typedef) as well as 
// C++ classes that might be used as base-classes for class declarations in the interface
// ie like above, this is how we get chbody stuff in i think

// %include "swig_test.h" // cannot do this because of vector library


void test();
// enum is weird, to change see https://stackoverflow.com/questions/16471213/wrapping-c-enum-in-a-python-module-with-swig
enum WaveMode { NONE, REGULAR }; // eventually add irregular waves mode
class HydroInputs {
public:
	WaveMode mode;
	HydroInputs();
	~HydroInputs() = default;
	double freq_index_des;
	double regular_wave_amplitude;
	double regular_wave_omega;
	double wave_omega_delta;
	std::vector<double> excitation_force_mag;
	std::vector<double> excitation_force_phase;
	HydroInputs(const HydroInputs& old);
	HydroInputs& operator = (const HydroInputs& rhs);
	void test2();
private:
};
class TestHydro {
public:
	bool printed = false;
	TestHydro();
	TestHydro(const std::vector<std::shared_ptr<ChBody>>& user_bodies, std::string h5_file_name, HydroInputs& users_hydro_inputs);
	TestHydro(const TestHydro& old) = delete;
	TestHydro operator = (const TestHydro& rhs) = delete;
	void WaveSetUp();
	std::vector<double> ComputeForceHydrostatics();
	std::vector<double> ComputeForceRadiationDampingConv();
	std::vector<double> ComputeForceExcitationRegularFreq();
	double GetRIRFval(int row, int col, int st);
	double coordinateFunc(int b, int i);
	bool convTrapz;
private:
	std::vector<std::shared_ptr<ChBody>> bodies;
	std::vector<H5FileInfo> file_info;
	std::vector<ForceFunc6d> force_per_body;
	double sumVelHistoryAndRIRF;
	HydroInputs hydro_inputs;
	std::vector<double> force_hydrostatic;
	std::vector<double> force_radiation_damping;
	std::vector<double> force_excitation_freq;
	std::vector<double> total_force;
	int num_bodies;
	std::vector<double> equilibrium;
	std::vector<double> cb_minus_cg;
	double rirf_timestep;
	double getVelHistoryAllBodies(int step, int c) const;
	double setVelHistory(double val, int step, int b_num, int index);
	std::vector<double> velocity_history; // use helper function to access vel_history elements correctly
	double prev_time;
	std::vector<double> rirf_time_vector; // (should be the same for each body?)
	int offset_rirf;
	std::shared_ptr<ChLoadContainer> my_loadcontainer;
	std::shared_ptr<ChLoadAddedMass> my_loadbodyinertia;
};