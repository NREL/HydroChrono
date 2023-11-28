Theory
======

Multibody dynamics theory in Project Chrono
-------------------------------------------

In this section, the fundamental multibody dynamics theory used by Chrono is described. For a multibody system comprising :math:`n_b` bodies, the generalized position of each body :math:`i` is represented by the vector :math:`r_i = [x_i, y_i, z_i]^T`, and its orientation is given by the Euler angles, :math:`e_i = [\psi_i, \theta_i, \phi_i]^T`. Hence, the generalized coordinates for the entire system can be expressed as:

.. math::
    :label: positions

    q = [r_1^T e_1^T \dots r_{n_b}^T e_{n_b}^T]^T \in \mathbb{R}^p,\ p = 6n_b

Where :math:`R` is the set of real numbers and :math:`p` is a variable representing the total number of generalized coordinates in the multibody system (i.e., each body has three variables for position and three for orientation).

In practice, quaternions (which have four variables) are typically used to describe orientation to avoid issues such as gimbal lock. If flexible bodies are present, deformation modes can be added to the position generalized coordinates.

Constrained mechanical systems have joints connecting bodies, which impose restrictions on the relative motion and result in constraints on the generalized coordinates. The kinematic constraints are formulated as algebraic expressions involving generalized coordinates:

.. math::
    :label: pos_constraints

    \Phi(q, t) = [\Phi_1(q, t) \dots \Phi_m(q, t)]^T = 0

where :math:`m` represents the total number of independent constraint equations that the generalized coordinates must satisfy throughout the simulation. For simplicity, only holonomic constraints (i.e., position-dependent constraints) are considered in this section.

Taking the time derivative of the equation yields the velocity kinematic constraint equation:

.. math::
    :label: vel_constraints

    \Phi_q(q, t)\dot{q} + \Phi_t(q, t) = 0

Here, :math:`\dot{q}` denotes generalized velocity, and the subscripts denote the partial derivatives: :math:`\Phi_q = \left[\frac{\partial \Phi_i}{\partial q_j}\right]` and :math:`\Phi_t = \left[\frac{\partial \Phi}{\partial t}\right]`. Differentiating the equation with respect to time results in the acceleration kinematic constraint equation:

.. math::
    :label: acc_constraints

    \Phi_q(q, t)\ddot{q} + (\Phi_q(q, t)\dot{q})q\dot{q} + 2\Phi_{qt}(q, t)\dot{q} + \Phi_{tt}(q, t) = 0

The time evolution of the mechanical system, governed by the Lagrange multiplier form of the constrained equations of motion, must satisfy the constraints at all times. The constrained equations of motion are given by:

.. math::
    :label: dae
    
    M(q)\ddot{q} + \Phi_q^T(q)\lambda = Q(\dot{q}, q, t)

where :math:`M(q) \in \mathbb{R}^{p \times p}` represents the generalized mass matrix, and :math:`Q(\dot{q}, q, t) \in \mathbb{R}^p` denotes the action force acting on the generalized coordinates :math:`q \in \mathbb{R}^p`. The reaction force :math:`\Phi_q^T(q)\lambda`, where :math:`\lambda \in \mathbb{R}^m` is the Lagrange multiplier associated with the kinematic constraints, results from the constraint equations. The combination of differential equations describing the system dynamics and algebraic equations describing the system's constraints leads to the system being referred to as a *differential-algebraic equation* (DAE). Because the constraint equations are described at the position level, these DAEs are referred to as "index-3."

The index-3 DAEs are neither linear nor ordinary differential. Hence, the Hilber-Hughes-Taylor (HHT) algorithm addresses these challenges by discretizing the equations of motion and enforcing position-level kinematic constraints. The core equation for the HHT integrator is:

.. math::
    :label: hht
    
    \frac{1}{1+\alpha}(M\ddot{q})_{n+1} + (\Phi^T_q \lambda - Q)_{n+1} - \frac{\alpha}{1+\alpha}(\Phi^T_q \lambda - Q)_n) = 0

The HHT algorithm iteratively solves the nonlinear equations for the unknowns using a Newton-like method, which eliminates ill-conditioning issues typically associated with integrating index-3 DAEs.

Hydrodynamics theory in HydroChrono
-----------------------------------

In order to integrate potential-flow-based hydrodynamics with Project Chrono, equation :eq:`dae` needs to be modified to incorporate added mass on its left side and the hydrodynamic force functions (e.g., hydrostatics, radiation damping, wave excitation) on its right side. Importantly, constraints, which might emerge from various subsystems like joints, power take-offs, and ropes, remain unchanged when including hydrodynamics. The time-integration of the system's equations of motion remains under Project Chrono's purview, without any need for alterations. Indeed, Chrono is solving the complete Cummins' equation.

The equation :eq:`dae` can be adapted to include the added mass and hydrodynamic forces as:

.. math::
    :label: dae_hydro

    (M(q) + A)\ddot{q} + \Phi_q^T(q)\lambda = Q(\dot{q}, q, t) + F_{hydro}(\dot{q}, q, t)

Where:
- :math:`A \in \mathbb{R}^{p \times p}` signifies the added mass matrix.
- :math:`F_{hydro}(\dot{q}, q, t) \in \mathbb{R}^p` represents the hydrodynamic force vector influencing the generalized coordinates :math:`q \in \mathbb{R}^p`.

The hydrodynamic force vector :math:`F_{hydro}(\dot{q}, q, t)` can be broken down into three primary components: hydrostatic force :math:`F_{hs}(q)`, radiation damping force :math:`F_{rad}(\dot{q}, q)`, and wave excitation force :math:`F_{exc}(t)`. Thus, we have:

.. math::
    :label: f_hydro
    
    F_{hydro}(\dot{q}, q, t) = F_{hs}(q, t) + F_{rad}(\dot{q}, t) + F_{exc}(t)

The subsequent sections provide an in-depth look into these components.

Hydrostatic force, :math:`F_{hs}(q, t)`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This force mirrors the restoring force due to buoyancy alongside alterations in the body's waterplane area. It can be derived from the product of the hydrostatic stiffness matrix :math:`K_h` and the displacement vector :math:`\Delta q` - the difference between the system's position :math:`q` and its equilibrium position :math:`q_0`:

.. math::
    :label: f_hs
    
    F_{hs}(q, t) = K_h \Delta q = K_h (q - q_0)

The hydrostatic stiffness matrix, :math:`K_h`, can be sourced by executing a numerical integration over the waterplane area of the floating body. Potential tools for this task include open-source boundary element method (BEM) tools like Capytaine :cite:`AnDi2019`, HAMS [jmse7030081]_, or the open-source mesh and hydrostatics package MeshMagick [#]_.

.. _[#]: https://github.com/LHEEA/meshmagick

Radiation damping force, :math:`F_{rad}(\dot{q}, t)`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This force stands as a representation of the energy dissipated when a floating body undergoes oscillation in water. It's procured through a convolution integral between the radiation impulse response function :math:`K_{rad}(t)` and the system's velocity timeline :math:`\dot{q}(\tau)`:

.. math::
    :label: f_rad
    
    F_{rad}(\dot{q}, t) = \int_{-\infty}^t K_{rad}(t - \tau) \dot{q}(\tau) d\tau

The :math:`K_{rad}(t)` function is derived by implementing the inverse continuous cosine transform (related to Fourier) on the frequency-domain radiation damping coefficients, :math:`B(\omega)`:

.. math::
    :label: k_rad
    
    K_{rad}(t) = \frac{2}{\pi} \int_0^\infty B(\omega) \cos(\omega t) d\omega

This transform allows the frequency domain coefficients, :math:`B(\omega)`, to be remapped into the time domain, thus producing the radiation impulse response function, :math:`K_{rad}(t)`. The :math:`B(\omega)` values can be sourced using open-source BEM software.

In HydroChrono, the force is computed through trapezoidal integration at the time values given by the RIRF time array relative to the current simulation time step. Linear interpolation is done for the velocity history if a given time value is between two values of the time series of the stored velocity history.
Wave excitation force, :math:`F_{exc}(t)`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following method to compute the wave excitation force involves convolution between the excitation Impulse Response Function (IRF) :math:`K_{exc}(t)` and the wave elevation time sequence :math:`\eta(t)`:

.. math::
    :label: f_ex

    F_{exc}(t) = \int_{-\infty}^{+\infty} K_{exc}(\tau) \eta(x, y, t-\tau) d\tau

By amalgamating these forces into the equation of motion, one can effectively model the behavior of a multibody oceanic system influenced by hydrodynamic forces.

In HydroChrono, the force is computed through trapezoidal integration by discretizing at the time values given by the excitation IRF time array relative to the current simulation time step. Linear interpolation is done for the free surface elevation if a given time value is between two values of the time series of the precomputed free surface elevation.


.. rubric:: Footnotes

.. [#] MeshMagick: https://github.com/LHEEA/meshmagick