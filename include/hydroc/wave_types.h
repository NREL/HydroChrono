/*********************************************************************
 * @file wave_types.h
 *
 * @brief header file for Wavebase and classes inheriting from WaveBase.
 *********************************************************************/
#pragma once
#include <hydroc/h5fileinfo.h>
#include <Eigen/Dense>

// TODO does M_PI need to be here?
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

// todo move this helper function somewhere else?
Eigen::VectorXd PiersonMoskowitzSpectrumHz(Eigen::VectorXd& f, double Hs, double Tp);
Eigen::VectorXd FreeSurfaceElevation(const Eigen::VectorXd& freqs_hz,
                                     const Eigen::VectorXd& spectral_densities,
                                     const Eigen::VectorXd& time_index,
                                     double water_depth,
                                     int seed = 1);
/**
 * @brief enum for type of wave.
 */
enum class WaveMode {
    /// @brief No waves
    noWaveCIC = 0,
    /// @brief Regular waves
    regular = 1,
    /// @brief Irregular waves
    irregular = 2
};

/**
 * @brief  pure virtual (interface) class for wave modes (regular, irregular, etc).
 */
class WaveBase {
  public:
    /**
     * @brief Override Initialize to set up any member variables in derrived class before timestepping occurs.
     */
    virtual void Initialize() = 0;

    /**
     * @brief Override to return the 6N dimensional force vector on hydro bodies.
     *
     * If force changes over time, put calculations
     *
     * @param t the current time to get the force for
     */
    virtual Eigen::VectorXd GetForceAtTime(double t) = 0;
    virtual WaveMode GetWaveMode()                   = 0;
};

/**
 * @brief class to instantiate WaveBase for no waves.
 */
class NoWave : public WaveBase {
  public:
    NoWave() { num_bodies = 1; }
    NoWave(unsigned int num_b) { num_bodies = num_b; }
    void Initialize() override {}

    /**
     * @brief fills a vector of 0s for force from NoWave.
     *
     * Overrides WaveBase function.
     *
     * @param t time, not used in NoWave case, here from WaveBase inheritance
     *
     * @return 6N dimensional force vector (Eigen::VectorXd) from waves in NoWave case.
     */
    Eigen::VectorXd GetForceAtTime(double t) override;
    WaveMode GetWaveMode() override { return mode; }

  private:
    unsigned int num_bodies;
    const WaveMode mode = WaveMode::noWaveCIC;
};

/**
 * @brief class to instantiate WaveBase for regular waves.
 */
class RegularWave : public WaveBase {
  public:
    /**
     * @brief default constructor for RegularWave.
     *
     * Defaults to 1 body if no argument given. Explicitly separate from constructor with num_b so a default \
     * constructor exists.
     */
    RegularWave();

    /**
     * @brief constructor for RegularWave in multibody cases.
     *
     * @param num_b number of bodies to apply hydro forces to (usually number of bodies in system)
     */
    RegularWave(unsigned int num_b);

    /**
     * @brief Initializes other member variables for timestep calculations later.
     *
     * Needs to be called before RegularWave can calculate any forces at any time. Resizes vectors, calculates \
     * constants, etc from h5 data. Assumes AddH5Data() has been called.
     */
    void Initialize() override;

    /**
     * @brief calculates the force from the regular wave at time t.
     *
     * Note if AdddH5Data and Initialize have not been called first, there will be issues.
     * TODO: add checks that these functions have been called in correct order.
     *
     * @param t time to calculate force at
     *
     * @return 6N dimensional Eigen::VectorXd for force from reg wave at time t, where N is number of bodies with \
     * hydroforces
     */
    Eigen::VectorXd GetForceAtTime(double t) override;

    /**
     * @brief gets wave mode.
     *
     * @return WaveMode enum for RegularWave is regular
     */
    WaveMode GetWaveMode() override { return mode; }

    // user input variables
    double regular_wave_amplitude;
    double regular_wave_omega;

    /**
     * @brief Initializes other member variables for timestep calculations later.
     *
     * Links the HydroData::RegularWaveInfo chunk to RegularWave for use in calculations.
     * Should be called before Initialize().
     *
     * @param reg_h5_data reference to chunk of h5 data needed for RegularWave calculations
     */
    void AddH5Data(std::vector<HydroData::RegularWaveInfo>& reg_h5_data);

  private:
    unsigned int num_bodies;
    const WaveMode mode = WaveMode::regular;
    std::vector<HydroData::RegularWaveInfo> wave_info;
    Eigen::VectorXd excitation_force_mag;
    Eigen::VectorXd excitation_force_phase;
    Eigen::VectorXd force;

    /**
     * @brief Finds omega_max and number of frequencies, then gets omega_max / num_freqs.
     *
     * Helper function for RegularWave calculations.
     *
     * @return omega step size
     */
    double GetOmegaDelta() const;

    /**
     * @brief gets interpolated excitation magnitude value.
     *
     * Helper function for RegWave calculations.
     *
     * @param b body to get matrix from
     * @param i row of matrix to look at
     * @param j column of matrix
     * @param freq_index_des TODO: what to call this
     *
     * @return excitation magnitudes for body b, row i, column j, frequency ix k
     */
    double GetExcitationMagInterp(int b, int i, int j, double freq_index_des) const;

    /**
     * @brief gets interpolated excitation phase value.
     *
     * Helper function for RegWave calculations.
     *
     * @param b body to get matrix from
     * @param i row of matrix to look at
     * @param j column of matrix
     * @param freq_index_des TODO: what to call this
     *
     * @return excitation phases for body b, row i, column j, frequency ix k
     */
    double GetExcitationPhaseInterp(int b, int i, int j, double freq_index_des) const;
};

/**
 * @brief class to instantiate WaveBase for irregular waves.
 */
class IrregularWave : public WaveBase {
  public:
    /**
     * @brief default constructor for IrregularWave.
     *
     * Defaults to 1 body if no argument given. Explicitly separate from constructor with num_b so a default \
     * constructor exists.
     */
    IrregularWave();

    /**
     * @brief constructor for IrregularWave in multibody case.
     *
     * @param num_b number of bodies to apply hydro forces to (usually number of bodies in system)
     */
    IrregularWave(unsigned int num_b);

    /**
     * @brief Initializes other member variables for timestep calculations later.
     *
     * Also resamples irf if timestep in h5 file does not match simulation timestep.
     * Needs to be called before IrregularWave can calculate any forces at any time. Resizes vectors, calculates \
     * constants, etc from h5 data. Assumes AddH5Data() has been called.
     * Calls CreateSpectrum() and CreateFreeSurfaceElevation()
     */
    void Initialize() override;

    /**
     * @brief computes 6N dimensional force from irregular wave on bodies in system.
     *
     * N is number of bodies with hydro forces (usually this is same as number of bodies in the system).
     *
     * @param t timestep to calculate the force at
     *
     * @return force vector (Eigen::VectorXd)
     */
    Eigen::VectorXd GetForceAtTime(double t) override;

    /**
     * @brief overloaded function from WaveBase to get the wave mode.
     *
     * @return WaveMode (always irregular for IrregularWave)
     */
    WaveMode GetWaveMode() override { return mode; }

    /**
     * @brief TODO
     *
     * @param start
     * @param end
     * @param num_steps
     *
     * @return eigen vector
     */
    Eigen::VectorXd SetSpectrumFrequencies(double start, double end, int num_steps);

    /**
     * @brief Makes an obj mesh for the current irregular wave.
     *
     * Computes points and triangles for free surface.
     * This is an optional function for users to call from their main simulation during set up if they want to have the
     * mesh exported and/or used in visualization.
     *
     * @param filename optionally set the file name to export mesh to. Must end in .obj for now. Defaults to
     * "fse_mesh.obj"
     */
    void SetUpWaveMesh(std::string filename = "fse_mesh.obj");

    /**
     * @brief Gives users access to the name of the mesh file generated by SetUpWaveMesh().
     *
     * @return file name from the wave's mesh .obj file.
     */
    std::string GetMeshFile();

    /**
     * @brief Gives users the wave mesh velocity for visualization purposes in simulation.
     *
     * wave mesh velocity is 1 in x direction
     *
     * @return wave mesh velocity as Eigen::Vector3d
     */
    Eigen::Vector3<double> GetWaveMeshVelocity();

    // user input // TODO add default values in case user doesn't initialize these?
    double wave_height;
    double wave_period;
    double simulation_duration;
    double simulation_dt;
    double ramp_duration;
    Eigen::VectorXd eta;  // public for mesh printing functions, TODO maybe make those friends, so this can be private?

    /**
     * @brief Initializes other member variables for timestep calculations later.
     *
     * Links the HydroData::IrregularWaveInfo chunk to IrregularWave for use in calculations.
     * Should be called before Initialize().
     *
     * @param irreg_h5_data reference to chunk of h5 data needed for IrregularWave calculations
     */
    void AddH5Data(std::vector<HydroData::IrregularWaveInfo>& irreg_h5_data, HydroData::SimulationParameters& sim_data);

  private:
    unsigned int num_bodies;
    const WaveMode mode = WaveMode::irregular;
    std::vector<HydroData::IrregularWaveInfo> wave_info;
    HydroData::SimulationParameters sim_data;
    std::vector<Eigen::MatrixXd> ex_irf_resampled;
    std::vector<Eigen::VectorXd> ex_irf_time_resampled;
    Eigen::VectorXd spectrum_frequencies;
    Eigen::VectorXd spectral_densities;
    std::string mesh_file_name;

    /**
     * @brief Get the excitation_irf_matrix from h5 file for specific body.
     *
     * @param b which body to get matrix for
     *
     * @return Eigen::MatrixXd irf matrix from h5 file for a specific body
     */
    Eigen::MatrixXd GetExcitationIRF(int b) const;

    /**
     * @brief resamples irf time vector from old time vector and new timestep.
     *
     * Creates a resized time vector that starts and ends at the same values as t_old and has timestep t_new.
     *
     * @param t_old reference to old time vector from h5 file, the first and last element are transfered to the first
     * and last element of t_new (return val)
     * @param dt_new the time step to use in resampled vector, typically the timestep of chrono simulation
     *
     * @return newly sized vector that looks like \
     * (t_old[0], t_old[0] + dt_new, t_old[0] + 2*dt_new, ... , t_old[t_old.size()-1])
     */
    Eigen::VectorXd ResampleTime(const Eigen::VectorXd& t_old, const double dt_new);

    /**
     * @brief Creates a resized values vector to interpolate values for a new timestep.
     *
     * @param t_old reference to old time vector from h5 file
     * @param vals_old the original values corresponding to the times in t_old from h5 file
     * @param t_new resampled times (return value from ResampleTime()) to use in interpolation
     *
     * @return matrix for interpolated vals_old over t_new
     */
    Eigen::MatrixXd ResampleVals(const Eigen::VectorXd& t_old, Eigen::MatrixXd& vals_old, const Eigen::VectorXd& t_new);

    /**
     * @brief Calculates the component of force from Convolution integral for specified body, dof, time.
     *
     * @param body which body currently calculating for
     * @param dof which degree of freedom to calculate force value for
     * @param time the time to compute force for
     *
     * @return value of force vector at t time in component corresponding to body and dof
     */
    double ExcitationConvolution(int body, int dof, double time);

    /**
     * @brief Creates wave spectrum for irrreg wave force.
     *
     * Currently only uses PiersonMoskowitzSpectrum, TODO make this adjustable for other spectra types.
     * Spectrum info saved to member variables for quick access in other functions.
     * Called from Initialize() function.
     */
    void CreateSpectrum();

    /**
     * @brief TODO describe
     *
     * Initializes eta and accounts for ramp_duration if set.
     * relevant info stored in member variables.
     */
    void CreateFreeSurfaceElevation();

    friend Eigen::VectorXd PiersonMoskowitzSpectrumHz(Eigen::VectorXd& f, double Hs, double Tp);
};

/**
 * @brief Opens obj file and writes points and triangles to file.
 *
 * @param points vector of tuples representing vertices of triangles
 * @param triangles list of faces each face is 3 indices of points
 * @param file_name name of obj file to write points and triangles to
 */
void WriteFreeSurfaceMeshObj(const std::vector<std::array<double, 3>>& points,
                             const std::vector<std::array<size_t, 3>>& triangles,
                             const std::string& file_name);

/**
 * @brief Makes triangles for mesh (visible from both sides).
 *
 * @param eta TODO
 *
 * @return vector of tuples representing triangles for mesh
 */
std::vector<std::array<size_t, 3>> CreateFreeSurfaceTriangles(size_t eta_size);

/**
 * @brief vertices for mesh.
 *
 * @param eta TODO
 * @param t_vec vector of times at each timestep in simulation, evenly spaced in most cases
 *
 * @return vector of tuples representing vertices for mesh
 */
std::vector<std::array<double, 3>> CreateFreeSurface3DPts(const Eigen::VectorXd& eta, const Eigen::VectorXd& t_vec);
