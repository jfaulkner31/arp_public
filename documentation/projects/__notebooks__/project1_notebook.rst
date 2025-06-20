.. _proj1_jupyter_notebook:

NRE 7203 - Advanced Reactor Physics - Homework 1
================================================

JONATHONS WORK FOR HOMEWORK 1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code::

    # import relevant packages
    import matplotlib.pyplot as plt
    import numpy as np
    import copy

    from pointsource_sphere import PointSourceInSphere
    from pointsource_sphere import PlotSigT


    plt.rcParams['font.size'] = 16
    %matplotlib inline
    plt.rcParams['figure.figsize'] = [6, 4] # Set default figure size

.. code::

    # Define some constants
    nMC = 100
    S0 = 5000
    R = 12.0
    nRegions_a = 10
    nRegions_b = 10
    nRegions_c = 30

    sigT_a = np.random.uniform(0.10, 0.25, size=nRegions_a)

    sigT_b = np.random.uniform(0.10, 0.25, size=nRegions_b)
    sigT_b[6] = 3.0 # manually assign to 7th pos in array (starting from 0 so really pos 6 in python world).

    sigT_c = np.random.uniform(0.10, 0.25, size=nRegions_c)



Calculation a
'''''''''''''

.. code::

    # Calculation a
    mc = PointSourceInSphere(nMC=nMC, S0=S0, R=R, sigT=sigT_a)
    mc.Solve(scheme='analytic')
    mc.Solve('ST')
    mc.Solve('DT')
    mc.Solve('wdt')

    # EXECUTION TIMES
    print("EXECUTION TIMES")
    print("ST = {:.5f}".format(mc.times['ST']))
    print("DT = {:.5f}".format(mc.times['DT']))
    print("WDT = {:.5f}".format(mc.times['WDT']))

    # LEAKAGES
    print("\nLEAKAGE")
    print("Analytic Leakage = {:.5f} %".format(100*mc.resAN['leakage']))
    print("MC/ST Leakage = {:.5f} %".format(100*mc.resST['leakage']), "pm {:.5f}".format(100*mc.resST['errleakage']))
    print("MC/DT Leakage = {:.5f} %".format(100*mc.resDT['leakage']), "pm {:.5f}".format(100*mc.resDT['errleakage']))
    print("MC/WDT Leakage = {:.5f} %".format(100*mc.resWDT['leakage']), "pm {:.5f}".format(100*mc.resWDT['errleakage']))


    # FIGURE OF MERITS
    print("\nFOMs")
    mc.PrintFOM('ST')
    mc.PrintFOM('DT')
    mc.PrintFOM('WDT')

    # FLUXES
    mc.PlotFluxes('all')

    # Deepcopy data:
    mc_a = copy.deepcopy(mc)

    # Save figure
    plt.savefig('Results/Figure_a.png')


.. parsed-literal::




    Now solving with method: analytic
    Now solving with method: ST
    Running generation 0
    Running generation 10
    Running generation 20
    Running generation 30
    Running generation 40
    Running generation 50
    Running generation 60
    Running generation 70
    Running generation 80
    Running generation 90
    Finished after 99 MC runs! Time was 3.873950958251953



    Now solving with method: DT
    Running generation 0
    Running generation 10
    Running generation 20
    Running generation 30
    Running generation 40
    Running generation 50
    Running generation 60
    Running generation 70
    Running generation 80
    Running generation 90
    Finished after 99 MC runs! Time was 3.5577402114868164



    Now solving with method: wdt
    Running generation 0
    Running generation 10
    Running generation 20
    Running generation 30
    Running generation 40
    Running generation 50
    Running generation 60
    Running generation 70
    Running generation 80
    Running generation 90
    Finished after 99 MC runs! Time was 4.57941746711731



    EXECUTION TIMES
    ST = 3.87395
    DT = 3.55774
    WDT = 4.57942

    LEAKAGE
    Analytic Leakage = 12.79313 %
    MC/ST Leakage = 12.81720 % pm 0.46496
    MC/DT Leakage = 12.86280 % pm 0.55292
    MC/WDT Leakage = 12.70944 % pm 0.39017

    FOMs
    FOM - ST: 0.00031194857915035066
    FOM - DT: 0.0002407959257579548
    FOM - WDT: 0.0003618540003301058



.. image:: ../NRE7203_HOMEWORK1_files/NRE7203_HOMEWORK1_5_1.png


Calculation b
'''''''''''''

.. code::

    mc = PointSourceInSphere(nMC=nMC, S0=S0, R=R, sigT=sigT_b)
    mc.Solve(scheme='analytic')
    mc.Solve('ST')
    mc.Solve('DT')
    mc.Solve('wdt')

    # EXECUTION TIMES
    print("EXECUTION TIMES")
    print("ST = {:.5f}".format(mc.times['ST']))
    print("DT = {:.5f}".format(mc.times['DT']))
    print("WDT = {:.5f}".format(mc.times['WDT']))

    # LEAKAGES
    print("Analytic Leakage = {:.5f} %".format(100*mc.resAN['leakage']))
    print("MC/ST Leakage = {:.5f} %".format(100*mc.resST['leakage']), "pm {:.5f}".format(100*mc.resST['errleakage']))
    print("MC/DT Leakage = {:.5f} %".format(100*mc.resDT['leakage']), "pm {:.5f}".format(100*mc.resDT['errleakage']))
    print("MC/WDT Leakage = {:.5f} %".format(100*mc.resWDT['leakage']), "pm {:.5f}".format(100*mc.resWDT['errleakage']))

    # FIGURE OF MERITS
    print("\nFOMs")
    mc.PrintFOM('ST')
    mc.PrintFOM('DT')
    mc.PrintFOM('WDT')

    # FLUXES
    mc.PlotFluxes('all')

    # Deepcopy data:
    mc_b = copy.deepcopy(mc)

    # Save figure
    plt.savefig('Results/Figure_b.png')



.. parsed-literal::




    Now solving with method: analytic
    Now solving with method: ST
    Running generation 0
    Running generation 10
    Running generation 20
    Running generation 30
    Running generation 40
    Running generation 50
    Running generation 60
    Running generation 70
    Running generation 80
    Running generation 90
    Finished after 99 MC runs! Time was 3.795607089996338



    Now solving with method: DT
    Running generation 0
    Running generation 10
    Running generation 20
    Running generation 30
    Running generation 40
    Running generation 50
    Running generation 60
    Running generation 70
    Running generation 80
    Running generation 90
    Finished after 99 MC runs! Time was 17.01483416557312



    Now solving with method: wdt
    Running generation 0
    Running generation 10
    Running generation 20
    Running generation 30
    Running generation 40
    Running generation 50
    Running generation 60
    Running generation 70
    Running generation 80
    Running generation 90
    Finished after 99 MC runs! Time was 35.381208419799805



    EXECUTION TIMES
    ST = 3.79561
    DT = 17.01483
    WDT = 35.38121
    Analytic Leakage = 2.90537 %
    MC/ST Leakage = 2.90860 % pm 0.21444
    MC/DT Leakage = 2.94320 % pm 0.22673
    MC/WDT Leakage = 2.89972 % pm 0.19785

    FOMs
    FOM - ST: 0.00035371562390919915
    FOM - DT: 8.012064679268595e-05
    FOM - WDT: 7.527196584226644e-05



.. image:: ../NRE7203_HOMEWORK1_files/NRE7203_HOMEWORK1_7_1.png


Calculation c
'''''''''''''

.. code::

    mc = PointSourceInSphere(nMC=nMC, S0=S0, R=R, sigT=sigT_c)
    mc.Solve(scheme='analytic')
    mc.Solve('ST')
    mc.Solve('DT')
    mc.Solve('wdt')

    # EXECUTION TIMES
    print("EXECUTION TIMES")
    print("ST = {:.5f}".format(mc.times['ST']))
    print("DT = {:.5f}".format(mc.times['DT']))
    print("WDT = {:.5f}".format(mc.times['WDT']))

    # LEAKAGES
    print("Analytic Leakage = {:.5f} %".format(100*mc.resAN['leakage']))
    print("MC/ST Leakage = {:.5f} %".format(100*mc.resST['leakage']), "pm {:.5f}".format(100*mc.resST['errleakage']))
    print("MC/DT Leakage = {:.5f} %".format(100*mc.resDT['leakage']), "pm {:.5f}".format(100*mc.resDT['errleakage']))
    print("MC/WDT Leakage = {:.5f} %".format(100*mc.resWDT['leakage']), "pm {:.5f}".format(100*mc.resWDT['errleakage']))

    # FIGURE OF MERITS
    print("\nFOMs")
    mc.PrintFOM('ST')
    mc.PrintFOM('DT')
    mc.PrintFOM('WDT')

    # FLUXES
    mc.PlotFluxes('all')

    # Deepcopy data:
    mc_c = copy.deepcopy(mc)

    # Save figure
    plt.savefig('Results/Figure_c.png')


.. parsed-literal::




    Now solving with method: analytic
    Now solving with method: ST
    Running generation 0
    Running generation 10
    Running generation 20
    Running generation 30
    Running generation 40
    Running generation 50
    Running generation 60
    Running generation 70
    Running generation 80
    Running generation 90
    Finished after 99 MC runs! Time was 6.552052974700928



    Now solving with method: DT
    Running generation 0
    Running generation 10
    Running generation 20
    Running generation 30
    Running generation 40
    Running generation 50
    Running generation 60
    Running generation 70
    Running generation 80
    Running generation 90
    Finished after 99 MC runs! Time was 3.69978404045105



    Now solving with method: wdt
    Running generation 0
    Running generation 10
    Running generation 20
    Running generation 30
    Running generation 40
    Running generation 50
    Running generation 60
    Running generation 70
    Running generation 80
    Running generation 90
    Finished after 99 MC runs! Time was 4.0094358921051025



    EXECUTION TIMES
    ST = 6.55205
    DT = 3.69978
    WDT = 4.00944
    Analytic Leakage = 10.48701 %
    MC/ST Leakage = 10.54000 % pm 0.41469
    MC/DT Leakage = 10.50660 % pm 0.42146
    MC/WDT Leakage = 10.44587 % pm 0.35399

    FOMs
    FOM - ST: 0.0002217139946331617
    FOM - DT: 0.00036997018345799365
    FOM - WDT: 0.000429480888408497



.. image:: ../NRE7203_HOMEWORK1_files/NRE7203_HOMEWORK1_9_1.png


PLOTTING SIGMA TOTAL
^^^^^^^^^^^^^^^^^^^^

.. code::

    PlotSigT(mc=mc_a, yLower=0.0, yUpper=0.3)



.. image:: ../NRE7203_HOMEWORK1_files/NRE7203_HOMEWORK1_11_0.png


.. code::

    PlotSigT(mc=mc_b, yLower=0.0, yUpper=3.1)



.. image:: ../NRE7203_HOMEWORK1_files/NRE7203_HOMEWORK1_12_0.png


.. code::

    PlotSigT(mc=mc_c, yLower=0.0, yUpper=0.3)



.. image:: ../NRE7203_HOMEWORK1_files/NRE7203_HOMEWORK1_13_0.png

