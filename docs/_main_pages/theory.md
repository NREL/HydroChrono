---
layout: page
title: Theory
permalink: /theory/
---

# Theory

## Multibody dynamics theory in [Project Chrono](https://projectchrono.org/)

In this section, the fundamental multibody dynamics theory used by Chrono is described. For a multibody system comprising \\( n_b \\) bodies, the generalized position of each body \\( i \\) is represented by the vector \\( r_i = [x_i, y_i, z_i]^T \\), and its orientation is given by the Euler angles, \\( e_i = [\psi_i, \theta_i, \phi_i]^T \\). Hence, the generalized coordinates for the entire system can be expressed as:

$$
q = [r_1^T e_1^T ... r_{n_b}^T e_{n_b}^T]^T ∈ ℝ^p, p = 6n_b
$$

Where \\( ℝ \\) is the set of real numbers and \\( p \\) is a variable representing the total number of generalized coordinates in the multibody system (i.e., each body has three variables for position and three for orientation).

In practice, quaternions (which have four variables) are typically used to describe orientation to avoid issues such as gimbal lock. If flexible bodies are present, deformation modes can be added to the position generalized coordinates.

Constrained mechanical systems have joints connecting bodies, which impose restrictions on the relative motion and result in constraints on the generalized coordinates. The kinematic constraints are formulated as algebraic expressions involving generalized coordinates:

$$
\Phi(q, t) = [\Phi_1(q, t) ... \Phi_m(q, t)]^T = 0
$$

where \\( m \\) represents the total number of independent constraint equations that the generalized coordinates must satisfy throughout the simulation. For simplicity, only holonomic constraints (i.e., position-dependent constraints) are considered in this section.

Taking the time derivative of the equation yields the velocity kinematic constraint equation:

$$
\Phi_q(q, t)\dot{q} + \Phi_t(q, t) = 0
$$

Here, \\( \dot{q} \\) denotes generalized velocity, and the subscripts denote the partial derivatives: \\( \Phi_q = [\frac{\partial \Phi_i}{\partial q_j}] \\) and \\( \Phi_t = [\frac{\partial \Phi}{\partial t}] \\). Differentiating the equation with respect to time results in the acceleration kinematic constraint equation:

$$
\Phi_q(q, t)\ddot{q} + (\Phi_q(q, t)\dot{q})q\dot{q} + 2\Phi_{qt}(q, t)\dot{q} + \Phi_{tt}(q, t) = 0
$$

The time evolution of the mechanical system, governed by the Lagrange multiplier form of the constrained equations of motion, must satisfy the constraints at all times. The constrained equations of motion are given by:

$$
M(q)\ddot{q} + \Phi_q^T(q)\lambda = Q(\dot{q}, q, t)
$$

where \\( M(q) ∈ ℝ^{p × p} \\) represents the generalized mass matrix, and \\( Q(\dot{q}, q, t) ∈ ℝ^p \\) denotes the action force acting on the generalized coordinates \\( q ∈ ℝ^p \\). The reaction force \\( \Phi_q^T(q)\lambda \\), where \\( \lambda ∈ ℝ^m \\) is the Lagrange multiplier associated with the kinematic constraints, results from the constraint equations. The combination of differential equations describing the system dynamics and algebraic equations describing the system's constraints leads to the system being referred to as a *differential-algebraic equation* (DAE). Because the constraint equations are described at the position level, these DAEs are referred to as "index-3."

The index-3 DAEs are neither linear nor ordinary differential. Hence, the Hilber-Hughes-Taylor (HHT) algorithm addresses these challenges by discretizing the equations of motion and enforcing position-level kinematic constraints. The core equation for the HHT integrator is:

$$
\frac{1}{1+\alpha}(M\ddot{q})_{n+1} + (\Phi^T_q \lambda - Q)_{n+1} - \frac{\alpha}{1+\alpha}(\Phi^T_q \lambda - Q)_n) = 0
$$

The HHT algorithm iteratively solves the nonlinear equations for the unknowns using a Newton-like method, which eliminates ill-conditioning issues typically associated with integrating index-3 DAEs.

## Hydrodynamics theory in HydroChrono

In order to integrate potential-flow-based hydrodynamics with Project Chrono, the main multibody dynamics equation needs to be modified to incorporate added mass on its left-hand side and the hydrodynamic force functions (e.g., hydrostatics, radiation damping, wave excitation) on its right-hand side. Importantly, constraints, which might emerge from various subsystems like joints, power take-offs, and ropes, remain unchanged when including hydrodynamics. The time-integration of the system's equations of motion is handled Project Chrono, without any need for alterations. Hence, Chrono is solving the complete Cummins' equation.

The main multibody dynamics equation of motion can be adapted to include the added mass and hydrodynamic forces as:

$$
(M(q) + A)\ddot{q} + \Phi_q^T(q)\lambda = Q(\dot{q}, q, t) + F_{hydro}(\dot{q}, q, t)
$$

Where:
- \\( A \in \mathbb{R}^{p \times p} \\) signifies the added mass matrix.
- \\( F_{hydro}(\dot{q}, q, t) \in \mathbb{R}^p \\) represents the hydrodynamic force vector influencing the generalized coordinates \\( q \in \mathbb{R}^p \\).

The hydrodynamic force vector \\( F_{hydro}(\dot{q}, q, t) \\) can be broken down into three primary components: hydrostatic force \\( F_{hs}(q) \\), radiation damping force \\( F_{rad}(\dot{q}, q) \\), and wave excitation force \\( F_{exc}(t) \\). Thus, we have:

$$
F_{hydro}(\dot{q}, q, t) = F_{hs}(q, t) + F_{rad}(\dot{q}, t) + F_{exc}(t)
$$

### Hydrostatic force, \\( F_{hs}(q, t) \\)

This force handles the restoring force due to buoyancy alongside alterations in the body's waterplane area. It can be derived from the product of the hydrostatic stiffness matrix \\( K_h \\) and the displacement vector \\( \Delta q \\) - the difference between the system's position \\( q \\) and its equilibrium position \\( q_0 \\):

$$
F_{hs}(q, t) = K_h \Delta q = K_h (q - q_0)
$$

The hydrostatic stiffness matrix, `K_h`, can be sourced by executing a numerical integration over the waterplane area of the floating body. Potential tools for this task include open-source boundary element method (BEM) tools like [Capytaine](https://capytaine.github.io/stable/), [HAMS](https://github.com/YingyiLiu/HAMS), or the open-source mesh and hydrostatics package [MeshMagick](https://lheea.github.io/meshmagick/).

### Radiation damping force, \\( F_{rad}(\dot{q}, t) \\)

This force stands as a representation of the energy dissipated when a floating body undergoes oscillation in water. It's procured through a convolution integral between the radiation impulse response function \\( K_{rad}(t) \\) and the system's velocity timeline \\( \dot{q}(\tau) \\):

$$
F_{rad}(\dot{q}, t) = \int_{-\infty}^t K_{rad}(t - \tau) \dot{q}(\tau) d\tau
$$

The \\( K_{rad}(t) \\) function is derived by implementing the inverse continuous cosine transform on the frequency-domain radiation damping coefficients, \\( B(\omega) \\):

$$
K_{rad}(t) = \frac{2}{\pi} \int_0^\infty B(\omega) \cos(\omega t) d\omega
$$

In HydroChrono, the force is computed through trapezoidal integration at the time values given by the RIRF time array relative to the current simulation time step. Linear interpolation is done for the velocity history if a given time value is between two values of the time series of the stored velocity history.

### Wave excitation force, \\( F_{exc}(t) \\)

The wave excitation force is defined by the convolution of the excitation Impulse Response Function (IRF) \\( K_{exc}(\tau) \\) with the free-surface elevation \\( \eta(t) \\):

$$
F_{exc}(t) = \int_{-\infty}^{+\infty} K_{exc}(\tau) \, \eta(x, y, t-\tau) \, d\tau
$$

By amalgamating these forces into the equation of motion, one can effectively model the behavior of a multibody oceanic system influenced by hydrodynamic forces.

#### Numerical evaluation

The continuous convolution integral is discretized using the quadrature points \\( \tau_j \\) and widths \\( \Delta\tau_j \\) provided by the excitation IRF time array from the BEM data:

$$
F_{exc,d}(t) = \sum_j K_{exc,d}(\tau_j) \, \eta(t - \tau_j) \, \Delta\tau_j
$$

where \\( d \\) denotes a particular degree of freedom.

HydroChrono supports two evaluation strategies, depending on how the sea state is defined.

**Frequency-domain transfer function (spectrum-based waves).** When the irregular sea state is represented as a linear superposition of \\( N_f \\) wave components,

$$
\eta(t) = \sum_{i=1}^{N_f} A_i \cos(\varphi_i - \omega_i \, t)
$$

the double summation over IRF points and frequency components can be factored by swapping the summation order and applying the angle-sum identity. This yields:

$$
F_{exc,d}(t) = \sum_{i=1}^{N_f} A_i \left[ C_{d,i} \cos(\varphi_i - \omega_i \, t) \;-\; S_{d,i} \sin(\varphi_i - \omega_i \, t) \right]
$$

where the *excitation transfer function* coefficients are pre-computed once from the IRF and the spectrum:

$$
C_{d,i} = \sum_j K_{exc,d}(\tau_j) \, \Delta\tau_j \, \cos(\omega_i \, \tau_j), \qquad
S_{d,i} = \sum_j K_{exc,d}(\tau_j) \, \Delta\tau_j \, \sin(\omega_i \, \tau_j)
$$

This reduces the per-timestep cost from \\( \mathcal{O}(N_{irf} \times N_f) \\) to \\( \mathcal{O}(N_f) \\) — the same cost as a single free-surface elevation evaluation — while producing results identical to the time-domain convolution to machine precision. Crucially, because the wave elevation is evaluated analytically, this formulation has no dependency on the simulation time step or duration, enabling the use of adaptive time integrators.

During the ramp-up period (when a cosine or linear ramp modifies \\( \eta \\)), the frequency-domain factoring does not apply. HydroChrono falls back to a batched matrix evaluation that computes \\( \eta \\) at all IRF quadrature points simultaneously via pre-computed trigonometric matrices, keeping the ramp startup efficient.

**Interpolation-based evaluation (imported elevation time series).** When the free-surface elevation is provided as an external time series (e.g., from a wave tank measurement or a coupled model), the convolution is evaluated directly by linearly interpolating the stored \\( \eta \\) data at each quadrature point \\( t - \tau_j \\). This path retains the traditional \\( \mathcal{O}(N_{irf}) \\) per-timestep cost.

## Links to related open-source software:

- [Project Chrono](https://projectchrono.org/)

- [Capytaine](https://capytaine.github.io/stable/)

- [HAMS](https://github.com/YingyiLiu/HAMS)

- [MeshMagick](https://github.com/LHEEA/meshmagick)

<p align="center">
  <img src="{{ site.baseurl }}/assets/img/wave_animation2.gif" alt="Wave Energy" width="80%" />
</p>