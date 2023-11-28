#ifndef HYDRO_FORCES_H
#define HYDRO_FORCES_H

/*********************************************************************
 * @file  hydro_forces.h
 *
 * @brief Header file of TestHydro main class and helper classes \
 * ComponentFunc and ForceFunc6d.
 *********************************************************************/

// TODO: clean up include statements

// Standard includes
#include <cstdio>
#include <filesystem>

// Chrono library includes
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

// Chrono FEA includes
#include <chrono/fea/ChMeshFileLoader.h>

// Hydroc library includes
#include <hydroc/h5fileinfo.h>
#include <hydroc/wave_types.h>

using namespace chrono;
using namespace chrono::fea;

class ForceFunc6d;
class TestHydro;

class ComponentFunc : public ChFunction {
  public:
    /**
     * @brief Default constructor.
     * @note Do not use this constructor.
     */
    ComponentFunc();

    /**
     * @brief Copy constructor.
     */
    ComponentFunc(const ComponentFunc& old);

    /**
     * @brief Construct a new ComponentFunc object.
     *
     * @param b Pointer to ForceFunc6d member object.
     * @param i Index for which component this ComponentFunc object refers to. Corresponds
     * to the force's degree of freedom, where:
     * (0,1,2,3,4,5) -> (surge, sway, heave, roll, pitch, yaw).
     */
    ComponentFunc(ForceFunc6d* b, int i);

    /**
     * @brief Clone the object.
     * @note This is a required override since ComponentFunc inherits from ChFunction.
     * @return A cloned ComponentFunc object.
     */
    virtual ComponentFunc* Clone() const override;

    /**
     * @brief Calculate the force value for a specific time.
     *
     * @param x Time from simulation.
     * @return Force on body in the specified degree of freedom at time x.
     */
    virtual double Get_y(double x) const override;

  private:
    ForceFunc6d* base_;  ///< Pointer to the full 6D force on the body.
    int index_;          ///< Index representing force degree of freedom on the body.
};

/**
 * @brief Organizes the functional (time-dependent) forces in each degree of freedom (6 total) for a body.
 */
class ForceFunc6d {
  public:
    /**
     * @brief Initializes an array of ComponentFunc objects and pointers to each force/torque.
     */
    ForceFunc6d();

    /**
     * @brief Initializes hydro force info from H5FileInfo and the ChBody this force will be applied to.
     *
     * @param object The body in the system to which this 6-dimensional force is being applied.
     * @param all_hydro_forces_user The TestHydro class where the total force on all bodies is calculated.
     */
    ForceFunc6d(std::shared_ptr<ChBody> object, TestHydro* all_hydro_forces_user);

    /**
     * @brief Copy constructor that ensures the force is only added to a body once.
     *
     * @note Avoid using the default copy constructor. ApplyForceAndTorqueToBody() should only ever be applied once.
     *
     * @param old ForceFunc6d object to copy from.
     */
    ForceFunc6d(const ForceFunc6d& old);

    /**
     * @brief Calculates the force on a given degree of freedom.
     *
     * @param i Index corresponding to the degree of freedom. Assumed to be 0-5 only.
     * @return Value of force in the i-th degree of freedom for the associated body.
     */
    double CoordinateFunc(int i);

  private:
    /**
     * @brief Initializes force components.
     */
    void SetForce();

    /**
     * @brief Initializes torque components.
     */
    void SetTorque();

    /**
     * @brief Adds this force to the body's list of applied forces.
     *
     * @warning Ensure this function is not called multiple times for the same force.
     */
    void ApplyForceAndTorqueToBody();

    std::shared_ptr<ChBody> body_;                  ///< Pointer to the body this force is applied to.
    int b_num_;                                     ///< Body's index in the system. Currently 1-indexed.
    ComponentFunc forces_[6];                       ///< Forces for each degree of freedom.
    std::shared_ptr<ComponentFunc> force_ptrs_[6];  ///< Pointers to the forces.
    std::shared_ptr<ChForce> chrono_force_;         ///< Chrono force for the body.
    std::shared_ptr<ChForce> chrono_torque_;        ///< Chrono torque for the body.
    TestHydro* all_hydro_forces_;                   ///< Pointer to TestHydro for calculations.
};

class ChLoadAddedMass;

// TODO: Rename TestHydro for clarity, perhaps to HydroForces?
// TODO: Split TestHydro class from its helper classes for clearer code structure.
class TestHydro {
  public:
    TestHydro() = delete;

    /**
     * @brief Main constructor for initializing the TestHydro class.
     *
     * Sets up vector of bodies, h5 file info, and hydro inputs. If no waves are given,
     * this constructor defaults to using NoWave.
     *
     * @param user_bodies List of pointers to bodies for the hydro forces.
     * @param h5_file_name Name of the h5 file where hydro data is stored.
     * @param waves WaveBase object. Defaults to NoWave if not provided.
     */
    TestHydro(std::vector<std::shared_ptr<ChBody>> user_bodies,
              std::string h5_file_name,
              std::shared_ptr<WaveBase> waves = std::make_shared<NoWave>());

    // Deleted copy constructor and assignment operator for safety.
    TestHydro(const TestHydro& old) = delete;
    TestHydro& operator=(const TestHydro& rhs) = delete;

    /**
     * @brief Adds waves class to force calculations depending on if regular or irregular waves.
     *
     * Also initializes h5 data for wave force class.
     *
     * @param waves the specific WaveBase class to add to the system
     */
    void AddWaves(std::shared_ptr<WaveBase> waves);

    /**
     * @brief Computes the Hydrostatic stiffness force plus buoyancy force for a 6N dimensional system.
     *
     * @return 6N dimensional force for 6 DOF and N bodies in system.
     */
    std::vector<double> ComputeForceHydrostatics();

    /**
     * @brief Computes the Radiation Damping force with convolution history for a 6N dimensional system.
     *
     * The discretization uses the time series of the the RIRF relative to the current step.
     * Linear interpolation is done on the velocity history if time_sim-time_rirf is between two values of the time
     * history. Trapezoidal integration is used to compute the force.
     *
     * Time history is automatically added in this function (so it should only be called once per time step), and
     * history that is older than the maximum RIRF time value is automatically removed.
     *
     * @return 6N dimensional force for 6 DOF and N bodies in system.
     */
    std::vector<double> ComputeForceRadiationDampingConv();

    /**
     * @brief Computes the 6N dimensional force from any waves applied to the system.
     * @return 6N dimensional force for 6 DOF and N bodies in system (already Eigen type).
     */
    Eigen::VectorXd ComputeForceWaves();

    /**
     * @brief Fetches the RIRF value from the h5 file based on the provided indices.
     *
     * @param row Index representing the body number and DOF index [0,...,5,...6N-1] for rows of RIRF.
     * @param col Column index in RIRF matrix [0,...,5,...6N-1].
     * @param st Index representing the timestep in RIRF, usually in the range [0,...1000].
     *
     * @return The value from the h5 file RIRF matrix.
     */
    double GetRIRFval(int row, int col, int st);

    /**
     * @brief Calculates or retrieves the total force on a specific body in a particular degree of freedom.
     *
     * If the total force for the body and DOF was computed for the current timestep, it's retrieved.
     * Otherwise, the function calculates it. Note: Body index is 1-based here due to its origin from ForceFunc6d.
     *
     * @param b Body index (1-based due to source from ForceFunc6d).
     * @param i Degree of Freedom (DOF) index, ranging from [0,...5].
     *
     * @return Component of the force vector for body 'b' and DOF 'i'.
     */
    double CoordinateFuncForBody(int b, int i);

  private:
    // Class properties related to the body and hydrodynamics
    std::vector<std::shared_ptr<ChBody>> bodies_;
    int num_bodies_;
    HydroData file_info_;
    std::vector<ForceFunc6d> force_per_body_;
    std::shared_ptr<WaveBase> user_waves_;

    // Force components vectors
    std::vector<double> force_hydrostatic_;
    std::vector<double> force_radiation_damping_;
    Eigen::VectorXd force_waves_;
    std::vector<double> total_force_;  // Saved force per timestep to reduce redundant calculations

    // Additional properties related to equilibrium and hydrodynamics
    std::vector<double> equilibrium_;
    std::vector<double> cb_minus_cg_;
    Eigen::VectorXd rirf_time_vector;  // Assumed consistent for each body
    Eigen::VectorXd rirf_width_vector;

    // Properties for velocity history management and time tracking
    std::vector<std::vector<std::vector<double>>> velocity_history_;
    std::vector<double> time_history_;
    double prev_time;

    // Added mass related properties
    std::shared_ptr<ChLoadContainer> my_loadcontainer;
    std::shared_ptr<ChLoadAddedMass> my_loadbodyinertia;

    /**
     * @brief Fetches the velocity history for a specific DOF, body, and timestep.
     *
     * @param step The timestep index, ranging from [0,1,...,1000].
     * @param c The combined index for the body and DOF, where the order is first by body then by DOF.
     *
     * @return Velocity history for the specified DOF and body at the given timestep.
     */
    double GetVelHistoryVal(int step, int c) const;

    /**
     * @brief Updates the velocity history for a given timestep, body, and DOF.
     *
     * @param val The new value to set.
     * @param step The timestep index, ranging from [0,1,...,1000].
     * @param b_num The body number (1-based), indicating which body's data to update.
     * @param index The DOF index, ranging from [0,1,...,5].
     */
    double SetVelHistory(double val, int step, int b_num, int index);
};

#endif