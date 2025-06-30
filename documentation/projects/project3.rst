.. _proj3:

Diffusion Coefficients and Critical Spectrum
---------------------------------------------

======================
Quick Scrolling
======================
* :ref:`Description <proj3_description>`
* :ref:`Methodology <proj3_methods>`
* :ref:`Results <proj3_results>`
* :ref:`Summary <proj3_summary>`
* :ref:`Jupyter Notebook <proj3_jupyter>`
* :ref:`Methods and Classes <proj3_classes>`


.. _proj3_description:

This project involves the calculation of diffusion coefficients through a variety of methods. The out-scattering, in-scattering, flux-limited, and Hydrogen-corrected methods are covered.
Additionally, the :math:`B_1`, :math:`P_1`, and CM methods for computing the critical spectrum (and then a diffusion coefficient) are discussed.

.. _proj3_methods:

-------------
In-Scatter
-------------

The in-scatter method comes from the P1 equations and computes the diffusion coefficient as:

		.. math::
			D_g = \frac{1}{3} \left( \Sigma_{t,g} - \frac{\sum_{g'} \Sigma_{s1,g'g} J_{g'}  }{J_g} \right)^{-1}

The current is generally not known however.
The implementation in this code is to compute the current using a fictitious buckling term to represent the flux gradient and solve the following P1 equation:

		.. math::
			\frac{1}{3} \frac{d\phi_g}{dx} + \Sigma_{t,g}J_g = \sum_{g'} \Sigma_{s1,g'g} J_{g'}

		.. math::
			\frac{1}{3} (B^2)^{0.5} + \Sigma_{t,g}J_g = \sum_{g'} \Sigma_{s1,g'g} J_{g'}

		.. math::
			\left( \underline{\underline{\Sigma_{s1,g'g}}} - \underline{\underline{\Sigma_{t,g}}} \right)\underline{J} = \frac{(B^2)^{0.5}}{3}[1 ... 1]^T

-------------
Out-Scatter
-------------
The out-scatter approach is much simpler. It can be written as:

		.. math::
			\Sigma_{tr,g} \approx \Sigma_{t,g} - \sum_{g'}\Sigma_{s1,g'g}

		.. math::
			D_g = \frac{1}{3\Sigma_{tr,g}}

-------------
Flux-Limited
-------------
The flux-limited approach is similar to the in-scatter approach except the neutron current, which is usually unknown, is replaced by the neutron flux.
The idea is to replace the neutron current by the known neutron flux spectra.

		.. math::
			\Sigma_{tr,g} \approx \Sigma_{t,g} - \frac{\sum_{g'}\Sigma_{s1,g'g} \phi_{g'}}{\phi_g}

		.. math::
			D_g = \frac{1}{3\Sigma_{tr,g}}

-------------------
Hydrogen-Corrected
-------------------
The Hydrogen-corrected approach considers most of the anisotropy to belong to the moderator (which is generally true for LWR's).
Two calculations are normally run - a medium of pure moderating material (pure Hydrogen) and a typical fuel assembly (FA).

The approach subtracts the "bad" out-scattering based component of the homogenized transport xs and then adds a "good" in-scatter based xs.

First, the transport correction ratio, :math:`\tau_H` for pure Hydrogen is obtained using a method (potentially the in-scatter approach) suitable for accurate calculations of transport cross sections
or transport correction factors.

Then, the pure Hydrogen microscopic cross section can be computed from the out-scatter transport xs and the Hydrogen number density in the pure Hydrogen medium (:math:`\sigma_{tr}^H = \Sigma_{tr}^H / N_{inf}^{H}`)




.. _proj3_results:

.. _proj3_summary:

.. _proj3_jupyter:

.. _proj3_classes:
