// module name is just a useful name for organizing files/namespaces/classes
// each module in a collection is created with separate swig invocations
// when modules need to share information (class derived from another in separate module) it's a bit complex 
// see import module
%module hydro_py_module 

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
%}

// However, #include statements are ignored unless the -includeall command line option has been supplied.  (out of %{...%} block)
// Now list ANSI C/C++ declarations
// this part is like a preprcessor....

// %inlcude over # includes each file once (dont need include guards)
// To include another file into a SWIG interface, use the %include directive like this:
// this makes more wrapper for pointer.i in this swig 
//%include "pointer.i"
%include "std_vector.i"
// Instantiate templates used by example
namespace std {
   %template(DoubleVector) vector<double>;
}

// %import "foo.i" on the other hand would import info from foo.i but not make more wrapper code
// Such information generally includes type declarations (e.g., typedef) as well as 
// C++ classes that might be used as base-classes for class declarations in the interface
// ie like above, this is how we get chbody stuff in i think

// %include "swig_test.h"


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
	HydroInputs(HydroInputs& old);
	HydroInputs& operator = (HydroInputs& rhs);
	void test2();
private:
};