/*********************************************************************
 * @file  hydro_forces.h
 *
 * @brief header file of TestHydro main class and helper classes \
 * ComponentFunc and ForceFunc6d.
 *********************************************************************/
// TODO: clean up include statements
#include <cstdio>
#include <filesystem>

#include <chrono/solver/ChIterativeSolverLS.h>
#include <chrono/solver/ChSolverPMINRES.h>
#include <chrono/timestepper/ChTimestepper.h>
#include <chrono/timestepper/ChTimestepperHHT.h>

#include <chrono/physics/ChBody.h>
#include <chrono/physics/ChBodyEasy.h>
#include <chrono/physics/ChForce.h>
#include <chrono/physics/ChLoad.h>
#include <chrono/physics/ChLoadContainer.h>
#include <chrono/physics/ChLoadsBody.h>
#include <chrono/physics/ChSystemNSC.h>
#include <chrono/physics/ChSystemSMC.h>

#include <chrono/fea/ChMeshFileLoader.h>

#include <hydroc/h5fileinfo.h>
#include <hydroc/wave_types.h>

using namespace chrono;
using namespace chrono::fea;

class ForceFunc6d;
class TestHydro;

class ComponentFunc : public ChFunction {
  public:
    /**
     * @brief default constructor, dont use.
     */
    ComponentFunc();
    // TODO delete default constructor for ComponentFunc? Users shouldnt be touching
    // ComponentFunc class anyway

    /**
     * @brief copy constructor works like default copy constructor.
     */
    ComponentFunc(const ComponentFunc& old);
    // TODO: remove copy constructor definition since it should be the default

    /**
     * @brief sets pointer to ForceFunc6d member object and index for which component
     * this ComponentFunc object refers to.
     *
     * @param b which body in the system the force is being applied to, 0 indexed?
     * @param i integer corresponding to the force's degree of freedom \
     * (0,1,2,3,4,5) -> (surge, sway, heave, roll, pitch, yaw)
     */
    ComponentFunc(ForceFunc6d* b, int i);

    /**
     * @brief Required override function since ComponentFunc inherits from ChFunction.
     */
    virtual ComponentFunc* Clone() const override;

    /**
     * @brief calculates the value of the index-th coordinate of the hydro forces at each time.
     *
     * This is a required override function since ComponentFunc inherits from ChFunction.
     *
     * @param x time from simulation
     *
     * @return force on body in index-th degree of freedom at time x
     */
    virtual double Get_y(double x) const override;

  private:
    ForceFunc6d* base;  // pointer to full 6d force on the body
    int index;          // which force degree of freedom this object represents on the body
};

// =============================================================================
// ForceFunc6d organizes the functional (time dependent) forces in each DoF (6 total) for a body
class ForceFunc6d {
  public:
    /**
     * @brief initializes array of ComponentFunc objects and pointers to each force/torque.
     */
    ForceFunc6d();

    /**
     * @brief calls default constructor and initializes hydro force info \
     * from H5FileInfo, also initializes ChBody that this force will be applied to.
     *
     * @param object which body in system this 6 dimensional force is being applied to
     * @param all_hydro_forces_user gets TestHydro class where total force on all bodies \
     * in the system is calculated so the components can be passed to ForceFunc6d to apply
     */
    ForceFunc6d(std::shared_ptr<ChBody> object, TestHydro* all_hydro_forces_user);

    /**
     * @brief copy constructor should check to see if this force has been added to this body yet
     * if not, it should add it, if so it shouldnt add the force a second time.
     *
     * Do not use default copy constructor for ForceFunc6d at this time! ApplyForceAndTorqueToBody()
     * should only ever be applied once. If your forces are doubled, tripled, or more
     * (depending on number of bodies) this is something to check (force only being applied once).
     */
    ForceFunc6d(const ForceFunc6d& old);

    /**
     * @brief calculates vector component of the force vector on the body.
     *
     * @param i index corresponding to the degree of freedom to calculate force on body \
     * assumed to be 0,1,2,3,4,5 only
     *
     * @return value of force in i-th degree of freedom for the body this force is applied to
     */
    double coordinateFunc(int i);

  private:
    /**
     * @brief Used to initialize components of force (external ChForce pointers).
     */
    void SetForce();

    /**
     * @brief used to initialize components of torque (external ChForce pointer with TORQUE flag set).
     */
    void SetTorque();

    // TODO: remove ApplyForceAndTorqueToBody() from ForceFunc6d, instead make the operations in this function
    // reachable by users in the main simulation loop. Like all forces, add to body in the demo
    // this will avoid a lot of issues we've had to overcome
    /**
     * @brief adds this force to the body's list of applied forces
     *
     * Warning: everytime this is called, a force is applied to the body so be careful not to \
     * duplicate forces on accident!!.
     * TODO: make this function less risky, there shouldn't be a chance to apply the same force
     * multiple times.
     */
    void ApplyForceAndTorqueToBody();

    std::shared_ptr<ChBody> body;  // pointer to body this 6d force is being applied to
    int b_num;                     // 1 indexed number representing which body in system TODO: make 0 indexed
    ComponentFunc forces[6];
    std::shared_ptr<ComponentFunc> force_ptrs[6];
    std::shared_ptr<ChForce> chrono_force;
    std::shared_ptr<ChForce> chrono_torque;
    TestHydro* all_hydro_forces;  // pointer up to TestHydro object so coordinateFunc() knows where to go for calcs
};

class ChLoadAddedMass;

// TODO give TestHydro a better name, perhaps HydroForces ?
// TODO split TestHydro class from its helper classes into new file above for clearer code
class TestHydro {
  public:
    bool printed = false;
    TestHydro()  = delete;
    /**
     * @brief main constructor for TestHydro class, sets up vector of bodies, h5 file info, and hydro inputs.
     *
     * Also initializes many persistent variables for force calculations from h5 file.
     *
     * @param user_bodies list of pointers to bodies in chrono system that the hydro forces are being applied to \
     * must be added to system in same order as provided here
     * @param h5_file_name name of h5 file where hydro data is stored
     * @param waves optional Wavebase parameter if using regular or irregular waves
     */
    TestHydro(std::vector<std::shared_ptr<ChBody>> user_bodies,
              std::string h5_file_name,
              std::shared_ptr<WaveBase> waves);

    /**
     * @brief alternate constructor for TestHydro class, sets up vector of bodies, h5 file info, and hydro inputs.
     *
     * Also initializes many persistent variables for force calculations from h5 file.
     * If no waves are given in main constructor, this constructor makes NoWave object and calls main constructor \
     * with NoWave object.
     * TODO: default NoWave option doesn't work for multibody sytems, add user_bodies->size() to NoWave constructor \
     * to fix?
     *
     * @param user_bodies list of pointers to bodies in chrono system that the hydro forces are being applied to \
     * must be added to system in same order as provided here
     * @param h5_file_name name of h5 file where hydro data is stored
     */
    TestHydro(std::vector<std::shared_ptr<ChBody>> user_bodies, std::string h5_file_name)
        : TestHydro(user_bodies, h5_file_name, std::static_pointer_cast<WaveBase>(std::make_shared<NoWave>())) {}

    // should only ever have 1 TestHydro object in a system to calc all the forces on all hydro bodies
    // don't try copying or moving this, many issues will arrise.
    TestHydro(const TestHydro& old) = delete;
    TestHydro operator=(const TestHydro& rhs) = delete;

    /**
     * @brief Adds waves class to force calculations depending on if regular or irregular waves.
     *
     * Also initializes h5 data for wave force class.
     *
     * @param waves the specific WaveBase class to add to the system
     */
    void AddWaves(std::shared_ptr<WaveBase> waves);

    void WaveSetUp();  // TODO remove function, no longer in use

    /**
     * @brief Computes the 6N dimensional Hydrostatic stiffness force plus buoyancy force.
     *
     * TODO: change return type to Eigen::VectorXd.
     * TODO: make private, users shouldn't call from main
     *
     * @return 6N dimensional force for 6 dof and N bodies in system
     */
    std::vector<double> ComputeForceHydrostatics();

    /**
     * @brief computes the 6N dimensional Radiation Damping force with convolution history.
     *
     * TODO: change return type to Eigen::VectorXd.
     * TODO: make private, users shouldn't call from main
     *
     * @return 6N dimensional force for 6 dof and N bodies in system
     */
    std::vector<double> ComputeForceRadiationDampingConv();

    /**
     * @brief computes the 6N dimensional force from any waves applied to system.
     *
     * TODO: make private, users shouldn't call from main
     *
     * @return 6N dimensional force for 6 dof and N bodies in system (already Eigen type)
     */
    Eigen::VectorXd ComputeForceWaves();

    /**
     * @brief looksup the rirf value from the h5 file info for correct body given the row, col, and step.
     *
     * TODO: make private function, since this is mostly a helper for ComputeForceRadiationDampingConv?
     *
     * @param row encodes the body number and dof index [0,...,5,...6N-1] for rows of RIRF
     * @param col col in RIRF matrix [0,...,5,...6N-1]
     * @param st which time step in rirf ranges usually [0,...1000]
     *
     * @return matrix val from h5 file RIRF val
     */
    double GetRIRFval(int row, int col, int st);

    /**
     * @brief Called from ForceFunc6d, this function finds the total force on body b in degree of freedom i.
     *
     * Calls all ComputeForce... type functions to get total force. If the total force has already been calculated \
     * this timestep, funciton instead returns the saved force from when it was calculated for this timestep.
     * b is 1 indexed here since it is from ForceFunc6d!!!!!.
     *
     * TODO: make private function and have ForceFunc6d be friends? users shouldn't be using this function from main
     *
     * @param b body number, it is 1 indexed here since it comes from ForcFunc6d TODO: make 0 indexed
     * @param i index or Degree of Freedom (dof) ranges [0,...5]
     *
     * @return component of force vector for body b (1 indexed) and degree of freedom i
     */
    double coordinateFunc(int b, int i);

    bool convTrapz;         // currently unused, see ComputeForceRadiationDampingConv()
    Eigen::VectorXd t_irf;  // TODO is this used in this class? if not, remove this variable

  private:
    std::vector<std::shared_ptr<ChBody>> bodies;
    int num_bodies;
    HydroData file_info;
    std::vector<ForceFunc6d> force_per_body;  // vector of ForceFunc6d on each body
    double sumVelHistoryAndRIRF;              // TODO what is this variable? is it used?
    std::shared_ptr<WaveBase> user_waves;     // applied wave object on bodies
    std::vector<double>
        force_hydrostatic;  // TODO do these need to be class level vectors, or can they me moved to compute functions?
    std::vector<double> force_radiation_damping;  // TODO do these need to be class level vectors, or can they me moved
                                                  // to compute functions?
    Eigen::VectorXd
        force_waves;  // TODO do these need to be class level vectors, or can they me moved to compute functions?
    std::vector<double> total_force;  // needs to be class level to save force each timestep (only calcs forces once,
                                      // then pulls from here)
    std::vector<double> equilibrium;
    std::vector<double> cb_minus_cg;
    double rirf_timestep;

    /**
     * @brief finds and returns the component of velocity history for given step and dof.
     *
     * helper function for ComputeForceRadiationDampingConv(), use to access vel history.
     *
     * @param step [0,1,...,1000] (timesteps from h5 file, one velocity per step
     * @param c column [0,..,num_bodies-1,...,numbodies*6-1] (in order of bodies, iterates over dof for each body...3 \
     * bodies c would be one of [0,1,...,17])
     *
     * @return velocity history for specified degree of freedom, body, and step (dof and b both specified in c)
     */
    double getVelHistoryVal(int step, int c) const;

    /**
     * @brief sets velocity history for step, b_num (body number) and index (dof) to the given val.
     *
     * helper function for ComputeForceRadiationDampingConv(), use to set/change vel history.
     *
     * @param val value to set the requested element to
     * @param step [0,1,...,1000] (0 indexed, up to the final timestep in h5 file)
     * @param b_num [1,2,...,total_bodies] (1 indexed!!!, use body number in h5 file) TODO make 0 indexed
     * @param index [0,1,2,3,4,5] (0 indexed, always 0-5 for force+torque vector indexing)
     */
    double setVelHistory(double val, int step, int b_num, int index);

    // double freq_index_des;
    // int freq_index_floor;
    // double freq_interp_val;
    std::vector<double> velocity_history;  // use helper function to access vel_history elements correctly
    double prev_time;                      // used to keep track of if force has been calculated this step or not
    Eigen::VectorXd rirf_time_vector;      // (should be the same for each body?)
    int offset_rirf;                       // used in circular nature of velocity history for convolution integral
    std::shared_ptr<ChLoadContainer> my_loadcontainer;    // stuff for added mass
    std::shared_ptr<ChLoadAddedMass> my_loadbodyinertia;  // stuff for added mass
};
