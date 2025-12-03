/*********************************************************************
 * @file  radiation_rirf_convolution.h
 * @brief RIRF-based radiation damping via velocity history convolution.
 *********************************************************************/

#ifndef HYDRO_RADIATION_RIRF_CONVOLUTION_H
#define HYDRO_RADIATION_RIRF_CONVOLUTION_H

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <functional>
#include <vector>

namespace hydrochrono::hydro {

/**
 * @brief RIRF-based radiation damping convolution.
 * 
 * Manages velocity history buffers and performs convolution of RIRF kernels
 * with body velocity history to compute radiation damping forces.
 */
class RadiationRirfConvolution {
public:
    /**
     * @brief Constructor.
     * 
     * @param num_bodies Number of bodies in the system
     * @param rirf_steps Number of time steps in RIRF
     * @param rirf_time_vector Time vector for RIRF
     * @param rirf_width_vector Width vector for RIRF (for trapezoidal integration)
     */
    RadiationRirfConvolution(
        int num_bodies,
        int rirf_steps,
        const Eigen::VectorXd& rirf_time_vector,
        const Eigen::VectorXd& rirf_width_vector);

    /**
     * @brief Record current velocities and update history buffers.
     * 
     * @param simulation_time Current simulation time
     * @param body_velocities Velocities for all bodies [body][dof] where dof is 6-DOF
     */
    void RecordVelocities(double simulation_time, 
                         const std::vector<std::vector<double>>& body_velocities);

    /**
     * @brief Compute radiation damping forces via RIRF convolution.
     * 
     * @param get_rirf_val Function to get RIRF value: (row, col, step) -> value
     *                      For Baseline mode: raw RIRF from H5
     *                      For TaperedDirect mode: processed RIRF tensor
     * @param processed_rirf Optional processed RIRF tensors (for TaperedDirect mode)
     *                       If provided, get_rirf_val is ignored and processed_rirf is used
     * 
     * @return Radiation damping force vector [6N] for N bodies
     */
    std::vector<double> ComputeForces(
        std::function<double(int, int, int)> get_rirf_val,
        const std::vector<Eigen::Tensor<double, 3>>* processed_rirf = nullptr);

private:
    int num_bodies_;
    int rirf_steps_;
    int total_dofs_;
    Eigen::VectorXd rirf_time_vector_;
    Eigen::VectorXd rirf_width_vector_;

    // Velocity history buffers (multibody: one per body)
    std::vector<std::vector<std::vector<double>>> velocity_history_;
    std::vector<double> time_history_;

    // Helper functions
    void PruneHistory(double history_min_time);
    static void InterpolateVelocity6D(
        const std::vector<std::vector<double>>& velocity_history_body,
        size_t newer_index,
        double query_time,
        double older_time,
        double newer_time,
        double out_velocity[6]);
    static bool AdvanceToBracket(
        const std::vector<double>& time_history,
        size_t& index,
        double query_time);
};

}  // namespace hydrochrono::hydro

#endif  // HYDRO_RADIATION_RIRF_CONVOLUTION_H

