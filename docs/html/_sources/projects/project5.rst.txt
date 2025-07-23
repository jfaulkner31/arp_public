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

Then, when calling certain functions, a simple if statement is added based on whether or not the symbolic representation
is chosen. The following is from :class:`analytic_nodal_expansion.AFEN2D._CornerPoint`:

.. code::

                if self.symbolic:
                    if idx > 0: # at idx = 0 we dont have symbolic, just floating number 1.0
                        val = xsi.subs(self.x, xval).subs(self.y, yval)
                    else:
                        val = 1.0
                else:
                    val = xsi(xval,yval,C,k)

Other similar if statements are made where necessary.

---------------------
2D power calculation
---------------------

A 2D power calculation was implemented in :ref:`Project 4's <proj4>` NEM module as
well as in the newly developed AFEN module. The power calculation basically
sums over the energy groups as follows to compute the power density:

		.. math::
			P''' = \sum_{g=1}^G \kappa \Sigma_{f,g} \phi_g

The power may be a 2D representation as well as a 1D representation.

.. _project5_results:

===========
Results
===========

---------------------
Analytical Integrals
---------------------

Results are first shown for the analytic integrals. A simple test function was written
as part of then :class:`analytic_nodal_expansion.AFEN2D` class - :class:`analytical.AFEN2D.verifyBasisFunctions`

The output from this function is shown below and also :ref:`HERE <proj5_notebook_results_basis_functions>`. All integrals and evaluations of the basis
functions are exactly the same between the symbolic and analytical representations.

.. parsed-literal::

    BASIS FUNCTIONS
    Group 0 | key0 | ANY = 1.0 | SYM = 1.00000000000000 | Diff. = 0
    Group 0 | key1 | ANY = 5.1558738300831655 | SYM = 5.15587383008317 | Diff. = 0
    Group 0 | key2 | ANY = 5.251955345558114 | SYM = 5.25195534555811 | Diff. = 0
    Group 0 | key3 | ANY = 1.4580732741460756 | SYM = 1.45807327414608 | Diff. = 0
    Group 0 | key4 | ANY = 1.7680434589622103 | SYM = 1.76804345896221 | Diff. = 0
    Group 0 | key5 | ANY = 2.3385772645467684 | SYM = 2.33857726454677 | Diff. = 0
    Group 0 | key6 | ANY = 3.441519739973118 | SYM = 3.44151973997312 | Diff. = 0
    Group 0 | key7 | ANY = 2.515313444492246 | SYM = 2.51531344449225 | Diff. = 0
    Group 0 | key8 | ANY = 3.701609950064031 | SYM = 3.70160995006403 | Diff. = 0
    Group 1 | key0 | ANY = 1.0 | SYM = 1.00000000000000 | Diff. = 0
    Group 1 | key1 | ANY = 4485.935524514339 | SYM = 4485.93552451434 | Diff. = 0
    Group 1 | key2 | ANY = 4485.9356359738085 | SYM = 4485.93563597381 | Diff. = 0
    Group 1 | key3 | ANY = 47.35470217398589 | SYM = 47.3547021739859 | Diff. = 0
    Group 1 | key4 | ANY = 47.36525961067778 | SYM = 47.3652596106778 | Diff. = 0
    Group 1 | key5 | ANY = 3889.625804027583 | SYM = 3889.62580402758 | Diff. = 0
    Group 1 | key6 | ANY = 3902.1147010857094 | SYM = 3902.11470108571 | Diff. = 0
    Group 1 | key7 | ANY = 3889.645789671481 | SYM = 3889.64578967148 | Diff. = 0
    Group 1 | key8 | ANY = 3902.1347508999575 | SYM = 3902.13475089996 | Diff. = 0
    INTEGRAL X
    Group 0 | key0 | ANY = 20 | SYM = 20.0000000000000 | Diff. = 0
    Group 0 | key1 | ANY = 44.839484098657664 | SYM = 44.8394840986577 | Diff. = 7.10542735760100e-15
    Group 0 | key2 | ANY = 44.01917141474316 | SYM = 44.0191714147432 | Diff. = 7.10542735760100e-15
    Group 0 | key3 | ANY = 29.16146548292151 | SYM = 29.1614654829215 | Diff. = 0
    Group 0 | key4 | ANY = 35.360869179244204 | SYM = 35.3608691792442 | Diff. = 0
    Group 0 | key5 | ANY = 30.370132891672107 | SYM = 30.3701328916721 | Diff. = 0
    Group 0 | key6 | ANY = 44.693589318954125 | SYM = 44.6935893189541 | Diff. = 0
    Group 0 | key7 | ANY = 28.23620350666294 | SYM = 28.2362035066629 | Diff. = -3.55271367880050e-15
    Group 0 | key8 | ANY = 41.55323547495101 | SYM = 41.5532354749510 | Diff. = -7.10542735760100e-15
    Group 1 | key0 | ANY = 20 | SYM = 20.0000000000000 | Diff. = 0
    Group 1 | key1 | ANY = 9857.195767841471 | SYM = 9857.19576784147 | Diff. = 0
    Group 1 | key2 | ANY = 9857.195522925382 | SYM = 9857.19552292538 | Diff. = 0
    Group 1 | key3 | ANY = 947.0940434797178 | SYM = 947.094043479718 | Diff. = 0
    Group 1 | key4 | ANY = 947.3051922135556 | SYM = 947.305192213556 | Diff. = 0
    Group 1 | key5 | ANY = 12087.191305922875 | SYM = 12087.1913059229 | Diff. = -1.81898940354586e-12
    Group 1 | key6 | ANY = 12126.001128653186 | SYM = 12126.0011286532 | Diff. = -1.81898940354586e-12
    Group 1 | key7 | ANY = 12087.129199933222 | SYM = 12087.1291999332 | Diff. = -1.81898940354586e-12
    Group 1 | key8 | ANY = 12125.93882325224 | SYM = 12125.9388232522 | Diff. = -1.81898940354586e-12
    INTEGRAL Y
    Group 0 | key0 | ANY = 10 | SYM = 10.0000000000000 | Diff. = 0
    Group 0 | key1 | ANY = 51.558738300831656 | SYM = 51.5587383008317 | Diff. = 0
    Group 0 | key2 | ANY = 52.519553455581146 | SYM = 52.5195534555811 | Diff. = 0
    Group 0 | key3 | ANY = 15.09497917397982 | SYM = 15.0949791739798 | Diff. = 1.77635683940025e-15
    Group 0 | key4 | ANY = 12.448554698022278 | SYM = 12.4485546980223 | Diff. = 0
    Group 0 | key5 | ANY = 41.55323547495101 | SYM = 41.5532354749510 | Diff. = -7.10542735760100e-15
    Group 0 | key6 | ANY = 28.236203506662935 | SYM = 28.2362035066629 | Diff. = 0
    Group 0 | key7 | ANY = 44.693589318954125 | SYM = 44.6935893189541 | Diff. = 0
    Group 0 | key8 | ANY = 30.370132891672103 | SYM = 30.3701328916721 | Diff. = 3.55271367880050e-15
    Group 1 | key0 | ANY = 10 | SYM = 10.0000000000000 | Diff. = 0
    Group 1 | key1 | ANY = 44859.355245143386 | SYM = 44859.3552451434 | Diff. = 0
    Group 1 | key2 | ANY = 44859.356359738085 | SYM = 44859.3563597381 | Diff. = 0
    Group 1 | key3 | ANY = 104.07831820701843 | SYM = 104.078318207018 | Diff. = -1.42108547152020e-14
    Group 1 | key4 | ANY = 104.05511976443623 | SYM = 104.055119764436 | Diff. = 0
    Group 1 | key5 | ANY = 12125.93882325224 | SYM = 12125.9388232522 | Diff. = -1.81898940354586e-12
    Group 1 | key6 | ANY = 12087.129199933224 | SYM = 12087.1291999332 | Diff. = -3.63797880709171e-12
    Group 1 | key7 | ANY = 12126.001128653188 | SYM = 12126.0011286532 | Diff. = -3.63797880709171e-12
    Group 1 | key8 | ANY = 12087.191305922875 | SYM = 12087.1913059229 | Diff. = -1.81898940354586e-12
    INTEGRAL XY
    Group 0 | key0 | ANY = 800 | SYM = 800.000000000000 | Diff. = 0
    Group 0 | key1 | ANY = 0.0 | SYM = 0 | Diff. = 0
    Group 0 | key2 | ANY = 1760.7668565897266 | SYM = 1760.76685658973 | Diff. = 2.27373675443232e-13
    Group 0 | key3 | ANY = 0.0 | SYM = 0 | Diff. = 0
    Group 0 | key4 | ANY = 995.8843758417822 | SYM = 995.884375841782 | Diff. = 1.13686837721616e-13
    Group 0 | key5 | ANY = 0.0 | SYM = 0 | Diff. = 0
    Group 0 | key6 | ANY = 0.0 | SYM = 0 | Diff. = 0
    Group 0 | key7 | ANY = 0.0 | SYM = 0 | Diff. = 0
    Group 0 | key8 | ANY = 1363.7063877369098 | SYM = 1363.70638773691 | Diff. = -2.27373675443232e-13
    Group 1 | key0 | ANY = 800 | SYM = 800.000000000000 | Diff. = 0
    Group 1 | key1 | ANY = 0.0 | SYM = 0 | Diff. = 0
    Group 1 | key2 | ANY = 394287.82091701525 | SYM = 394287.820917015 | Diff. = 0
    Group 1 | key3 | ANY = 0.0 | SYM = 0 | Diff. = 0
    Group 1 | key4 | ANY = 8324.409581154898 | SYM = 8324.40958115490 | Diff. = 0
    Group 1 | key5 | ANY = 0.0 | SYM = 0 | Diff. = 0
    Group 1 | key6 | ANY = 0.0 | SYM = 0 | Diff. = 0
    Group 1 | key7 | ANY = 0.0 | SYM = 0 | Diff. = 0
    Group 1 | key8 | ANY = 150244.47045224538 | SYM = 150244.470452245 | Diff. = -2.91038304567337e-11

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
