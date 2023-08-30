###############################################################
Oscillating Surge Wave Energy Converter (OSWEC) - Verification
###############################################################

Overview
========

The Oscillating Surge Wave Energy Converter (OSWEC) is the subject of focused verification in this section, using the example model provided in the WEC-Sim package :cite:`WEC-SimReference`. The verification process serves as a comparative measure against WEC-Sim results for assessing the hydrodynamic performance and numerical reliability of the OSWEC model under study.

.. figure:: images/oswec_mesh_and_viz.png
   :scale: 25%
   :alt: Image

   Visualization of OSWEC. Left\: hydrodynamic mesh for coefficient computation. Right\: OSWEC system in the Chrono GUI.

Model Parameters
================

The OSWEC model's primary physical parameters—such as hinge location, center of gravity, and flap mass—are vital for a comprehensive hydrodynamic analysis.

.. raw:: html

    <style>
      .center-text {
        text-align: center;
      }
    </style>

    <table border="1" style="width:80%;">
      <thead>
        <tr>
          <th>Name</th>
          <th class="center-text">Symbol</th>
          <th class="center-text">Value</th>
          <th class="center-text">Units</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <td>Hinge Location</td>
          <td class="center-text">\[ \text{Hinge} \]</td>
          <td class="center-text">\[ \begin{bmatrix} 0.0 \\ 0.0 \\ -8.9 \end{bmatrix} \]</td>
          <td class="center-text">\[ \text{m} \]</td>
        </tr>
        <tr>
          <td>Flap Center of Gravity</td>
          <td class="center-text">\[ \text{CoG} \]</td>
          <td class="center-text">\[ \begin{bmatrix} 0.0 \\ 0.0 \\ -3.9 \end{bmatrix} \]</td>
          <td class="center-text">\[ \text{m} \]</td>
        </tr>
        <!-- Additional rows can be added similarly -->
      </tbody>
    </table>

Results
=======

The OSWEC model's numerical simulations undergo a verification process, particularly focusing on pitch decay tests and Response Amplitude Operators (RAOs).

.. figure:: images/oswec_decay_verification.png
   :scale: 25%
   :alt: Image

   OSWEC 10-degree pitch decay test, displaying a high degree of correlation between HydroChrono and WEC-Sim results.

.. figure:: images/oswec_rao_verification.png
   :scale: 25%
   :alt: Image

   RAOs from WEC-Sim and HydroChrono for OSWEC in regular waves, revealing minor discrepancies at the resonance frequency.

Conclusions
===========

The verification process substantiates the numerical implementation of the hydrodynamic forces for rotational degrees of freedom in the OSWEC model. Limitations remain, particularly in understanding hydrodynamic interactions in multibody systems, warranting further investigation.

References
==========

For extended insights into the verification process and the intricacies of the OSWEC model, consult :cite:`YourReference1` and :cite:`YourReference2`.
