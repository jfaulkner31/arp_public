.. _proj5:

Project 5: Analytic Flux Nodal Expansion Method
-----------------------------------------------

======================
Quick Scrolling
======================
* :ref:`Description <_project5_description>`
* :ref:`Methodology <project5_methods>`
* :ref:`Results <project5_results>`
* :ref:`Summary <project5_summary>`
* :ref:`Jupyter Notebook <project5_jupyter>`
* :ref:`Methods and Classes <project5_classes>`


.. _project5_description:

===========
Description
===========

This project involves the implementation of the
Analytic Flux Nodal Expansion Method (AFEN). This method typically takes boundary
conditions in the form of side fluxes or corner fluxes
from a reference (e.g. Serpent) or a NEM tool. AFEN provides a 2D representation
of the homogenous flux upon solving which can be used directly to compute the
2D homogenous power distribution. Normally, the 2D power is also paired with
form factors to better reconstruct heterogenous fluxes; this last part, however,
is not done in this project. (See Project 6)

.. _project5_methods:

===========
Methodology
===========

This project consists of three parts:

1. Implementation of an analytical method for computing integrals in place of Python's symbolic toolbox.
2. Calculation of the homogenous power density for each fuel assembly in the colorset.
3. Calculation of 1D averaged quantities with comparison to Serpent.

The analytic integrals and the power calculation are now discussed.

---------------------
Analytical Integrals
---------------------
The Python symbolic toolbox was previously used to compute certain integrals used
in the AFEN basis functions. To provide a massive boost to speedup, analytic functions
were pre-evaluated and input into the code in place of the symbolic toolbox.

The implementation was made so that the same data structures used previously may be
used again for the analytical integrals with only minimal changes to the code.

The following is a short code example from :class:`analytic_nodal_expansion.AFEN2D._IntegrateY`

Essentially, a flag "symbolic" has been added to the :class:`analytic_nodal_expansion.AFEN2D` class.
Based on user input, the symbolic or analytic expression is chosen. Python's lambda functions are used
to store the necessary functions in the same data structure previously used for the symbolic representation.


.. code::

    def _IntegrateY(self):
        """y integration over all the expansion functions"""

        xsi_int = {}
        if self.symbolic:
            for ig in range(self.ng):
                xsi_int[ig] = {}
                for idx, xsi in self.xsi[ig].items():
                    xsi_int[ig][idx] = sym.integrate(xsi, self.y)
        else:
            for ig in range(self.ng):
                lam = self.eigenvals[ig]
                xsi_int[ig] = {}
                if lam > 0:
                    xsi_int[ig][0] = lambda x,y,C,k: y
                    xsi_int[ig][1] = lambda x,y,C,k: np.sinh(k*x) * y
                    xsi_int[ig][2] = lambda x,y,C,k: np.cosh(k*x) * y
                    xsi_int[ig][3] = lambda x,y,C,k: np.cosh(k*y) / k
                    xsi_int[ig][4] = lambda x,y,C,k: np.sinh(k*y) / k
                    xsi_int[ig][5] = lambda x,y,C,k: np.sinh(k*C*x)*np.cosh(k*C*y) / k / C
                    xsi_int[ig][6] = lambda x,y,C,k: np.sinh(k*C*x)*np.sinh(k*C*y) / k / C
                    xsi_int[ig][7] = lambda x,y,C,k: np.cosh(k*C*x)*np.cosh(k*C*y) / k / C
                    xsi_int[ig][8] = lambda x,y,C,k: np.cosh(k*C*x)*np.sinh(k*C*y) / k / C
                else:
                    xsi_int[ig][0] = lambda x,y,C,k: y
                    xsi_int[ig][1] = lambda x,y,C,k: np.sin(k*x) * y
                    xsi_int[ig][2] = lambda x,y,C,k: np.cos(k*x) * y
                    xsi_int[ig][3] = lambda x,y,C,k: -np.cos(k*y) / k
                    xsi_int[ig][4] = lambda x,y,C,k: np.sin(k*y) / k
                    xsi_int[ig][5] = lambda x,y,C,k: np.sin(k*C*x)*-np.cos(k*C*y) / k / C
                    xsi_int[ig][6] = lambda x,y,C,k: np.sin(k*C*x)*np.sin(k*C*y) / k / C
                    xsi_int[ig][7] = lambda x,y,C,k: np.cos(k*C*x)*-np.cos(k*C*y) / k / C
                    xsi_int[ig][8] = lambda x,y,C,k: np.cos(k*C*x)*np.sin(k*C*y) / k / C

        self.xsiIntY = xsi_int

---------------------
2D power calculation
---------------------



---------------------
Analytical Integrals
---------------------

.. _project5_results:

===========
Results
===========

.. _project5_summary:

=================================
Summary
=================================

.. _project5_jupyter:

Jupyter Notebook
---------------------------------------------

.. _project5_classes:

Classes and Objects
---------------------------------------------
