.. _project6_notebook:

Final Project
~~~~~~~~~~~~~

.. code:: ipython3

    # import relevant packages
    import numpy as np
    from scipy.linalg import solve
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    import serpentTools
    from serpentTools.settings import rc

    from diffusion_coeffs import *

    import NEM

    import AFEN

    import time

    # Default values
    # rcParams['figure.dpi'] = 300
    plt.rcParams['figure.figsize'] = [6, 4] # Set default figure size

    # rc serpent version
    rc["serpentVersion"] = "2.2.1"

Question 1
~~~~~~~~~~

.. code:: ipython3

    # Question 1
    # part 1: condense from 70->2g for each FA in the FA onlly assemblies
    resFile_455 = "./serpent/Fuel_only/Fuel_455_70g_900K_res.m"
    resFile_260 = "./serpent/Fuel_only/Fuel_260_70g_res.m"

    res_455 = serpentTools.read(resFile_455)
    res_260 = serpentTools.read(resFile_260)


    #########
    # 4.55% #
    #########
    print("4.55% Below:")
    ng = 70
    univ0 = res_455.getUniv('0', timeDays=0)
    flx = univ0.infExp['infFlx']
    sigT = univ0.infExp['infTot']
    SP1 = univ0.infExp['infSp1']
    cmmTransp_FA = univ0.gc['cmmTranspxs']      # represents in-scatter
    infTransp = univ0.infExp['infTranspxs']  # out-scatter
    energy = univ0.groups * 1E+06
    SP1=SP1.reshape((ng,ng)).transpose()

    # Inscattering
    transportxs, tau_in_FA, Jg = InScatter(ng=70, sigS1=SP1, sigT=sigT, flx=flx, B2=0.0001) # get inscatter
    D = Condense2gr(xs=1/3/transportxs, flx=flx,energy=energy)
    print("Inscatter Diffusion Coeff. =", D)

    # Outscatter
    transportxs, tau_out_FA = OutScatter(ng=70, sigS1=SP1, sigT=sigT)
    D = Condense2gr(xs=1/3/transportxs, flx=flx,energy=energy)
    print("Outscatter Diffusion Coeff. =", D)

    # Flx Limited
    transportxs, tau_flx_FA = FluxLimited(ng=70, sigS1=SP1, sigT=sigT, flx=flx)
    D = Condense2gr(xs=1/3/transportxs, flx=flx,energy=energy)
    print("Flux Limited Diffusion Coeff. =", D)

    # Raw Condensation
    raw_condensed_diff = Condense2gr(1/3/cmmTransp_FA, flx, energy)
    raw_condensed_inv_diff = 1/3/Condense2gr(cmmTransp_FA, flx, energy)
    print("condensed CMM (1/3/tr method):",raw_condensed_diff)
    print("condensed CMM (tr method):",raw_condensed_inv_diff)

    # Reference


    #########
    # 2.60% #
    #########
    print("2.60% Below:")
    ng = 70
    univ0 = res_260.getUniv('0', timeDays=0)
    flx = univ0.infExp['infFlx']
    sigT = univ0.infExp['infTot']
    SP1 = univ0.infExp['infSp1']
    cmmTransp_FA = univ0.gc['cmmTranspxs']      # represents in-scatter
    infTransp = univ0.infExp['infTranspxs']  # out-scatter
    energy = univ0.groups * 1E+06
    SP1=SP1.reshape((ng,ng)).transpose()

    # Inscattering
    transportxs, tau_in_FA, Jg = InScatter(ng=70, sigS1=SP1, sigT=sigT, flx=flx, B2=0.0001) # get inscatter
    D = Condense2gr(xs=1/3/transportxs, flx=flx,energy=energy)
    print("Inscatter Diffusion Coeff. =", D)

    # Outscatter
    transportxs, tau_out_FA = OutScatter(ng=70, sigS1=SP1, sigT=sigT)
    D = Condense2gr(xs=1/3/transportxs, flx=flx,energy=energy)
    print("Outscatter Diffusion Coeff. =", D)

    # Flx Limited
    transportxs, tau_flx_FA = FluxLimited(ng=70, sigS1=SP1, sigT=sigT, flx=flx)
    D = Condense2gr(xs=1/3/transportxs, flx=flx,energy=energy)
    print("Flux Limited Diffusion Coeff. =", D)

    # Raw Condensation
    raw_condensed_diff = Condense2gr(1/3/cmmTransp_FA, flx, energy)
    raw_condensed_inv_diff = 1/3/Condense2gr(cmmTransp_FA, flx, energy)
    print("condensed CMM (1/3/tr method):",raw_condensed_diff)
    print("condensed CMM (tr method):",raw_condensed_inv_diff)





.. parsed-literal::

    4.55% Below:
    Inscatter Diffusion Coeff. = [1.36809841 0.39875562]
    Outscatter Diffusion Coeff. = [1.47303773 0.3621981 ]
    Flux Limited Diffusion Coeff. = [1.3210531  0.38653842]
    condensed CMM (1/3/tr method): [1.39156005 0.40622894]
    condensed CMM (tr method): [1.12251807 0.34999056]
    2.60% Below:
    Inscatter Diffusion Coeff. = [1.36078691 0.39423531]
    Outscatter Diffusion Coeff. = [1.46493777 0.3693202 ]
    Flux Limited Diffusion Coeff. = [1.31462087 0.38610637]
    condensed CMM (1/3/tr method): [1.38438031 0.40681814]
    condensed CMM (tr method): [1.11932516 0.35881331]


Question 2
~~~~~~~~~~

.. code:: ipython3

    ### interpolate the 4.55% from 600 + 1200 -> 900 K
    resFile_600 = "./serpent/Fuel_only/Fuel_455_2g_600K_res.m"
    resFile_1200 = "./serpent/Fuel_only/Fuel_455_2g_1200K_res.m"
    resFile_900 = "./serpent/Fuel_only/Fuel_455_2g_900K_res.m"

    res_600 = serpentTools.read(resFile_600)
    res_1200 = serpentTools.read(resFile_1200)
    res_900 = serpentTools.read(resFile_1200)

    ng = 2

    # 600
    univ0 = res_600.getUniv('0', timeDays=0)
    nsf_600 = univ0.infExp['infNsf']
    cmmTransp_FA_600 = univ0.gc['cmmTranspxs']      # represents in-scatter
    D_600 = 1/3/cmmTransp_FA_600
    SP0_600 = univ0.infExp['infSp0']
    SP0_600=SP0_600.reshape((ng,ng)).transpose()
    redabs_600 = univ0.infExp['infRabsxs']

    # 1200
    univ0 = res_1200.getUniv('0', timeDays=0)
    nsf_1200 = univ0.infExp['infNsf']
    cmmTransp_FA_1200 = univ0.gc['cmmTranspxs']      # represents in-scatter
    D_1200 = 1/3/cmmTransp_FA_1200
    SP0_1200 = univ0.infExp['infSp0']
    SP0_1200=SP0_1200.reshape((ng,ng)).transpose()
    redabs_1200 = univ0.infExp['infRabsxs']

    # 900
    univ0 = res_900.getUniv('0', timeDays=0)
    nsf_900 = univ0.infExp['infNsf']
    cmmTransp_FA_900 = univ0.gc['cmmTranspxs']      # represents in-scatter
    D_900 = 1/3/cmmTransp_FA_900
    SP0_900 = univ0.infExp['infSp0']
    SP0_900=SP0_900.reshape((ng,ng)).transpose()
    redabs_900 = univ0.infExp['infRabsxs']

    print("interpolatead D:", interpolateXS(x1=600, x2=1200, y1=D_600, y2=D_1200, x=900, interp_type='loglog'))
    print("interpolatead NSF:", interpolateXS(x1=600, x2=1200, y1=nsf_600, y2=nsf_1200, x=900, interp_type='loglog'))
    print("interpolatead SP0:", interpolateXS(x1=600, x2=1200, y1=SP0_600, y2=SP0_1200, x=900, interp_type='loglog'))
    print("interpolatead redabs:", interpolateXS(x1=600, x2=1200, y1=redabs_600, y2=redabs_1200, x=900, interp_type='loglog'))

    print("ref D", D_900)
    print("ref nsf", nsf_900)
    print("ref redabs", redabs_900)
    print("ref sp0", SP0_900)


.. parsed-literal::

    interpolatead D: [1.39151829 0.40634751]
    interpolatead NSF: [0.00863624 0.17273652]
    interpolatead SP0: [[0.52972432 0.00238746]
     [0.01711861 1.27196665]]
    interpolatead redabs: [0.01100566 0.12497814]
    ref D [1.39175939 0.40653133]
    ref nsf [0.00863375 0.172371  ]
    ref redabs [0.0111302 0.124775 ]
    ref sp0 [[0.529787   0.00256537]
     [0.0170528  1.27156   ]]


Question 3
~~~~~~~~~~

Part 1 - report corner DFs
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    resFile_455 = "./serpent/Fuel_only/Fuel_455_2g_900K_res.m"
    resFile_260 = "./serpent/Fuel_only/Fuel_260_2g_res.m"
    xs_455, bc_455 = NEM.GetSerpentRes(resFile=resFile_455, univ='0', timeDays=0)
    xs_260, bc_260 = NEM.GetSerpentRes(resFile=resFile_260, univ='0', timeDays=0)


.. code:: ipython3

    bc_455




.. parsed-literal::

    {'area': array([21.42, 21.42, 21.42, 21.42]),
     'flux': array([2.37961e+14, 3.19833e+13]),
     'wJnet': array([0., 0.]),
     'sJnet': array([0., 0.]),
     'eJnet': array([-0., -0.]),
     'nJnet': array([-0., -0.]),
     'wFlux': array([2.41842e+14, 3.37817e+13]),
     'sFlux': array([2.41769e+14, 3.37144e+13]),
     'eFlux': array([2.41827e+14, 3.37649e+13]),
     'nFlux': array([2.41898e+14, 3.37753e+13]),
     'wDF': array([1.01631, 1.05623]),
     'sDF': array([1.016  , 1.05413]),
     'eDF': array([1.01624, 1.0557 ]),
     'nDF': array([1.01654, 1.05603]),
     'nwDF': array([1.02128, 1.02203]),
     'neDF': array([1.02252, 1.02245]),
     'seDF': array([1.02076, 1.02134]),
     'swDF': array([1.02235, 1.02145])}



.. code:: ipython3

    bc_260




.. parsed-literal::

    {'area': array([21.42, 21.42, 21.42, 21.42]),
     'flux': array([2.38455e+14, 5.33629e+13]),
     'wJnet': array([0., 0.]),
     'sJnet': array([0., 0.]),
     'eJnet': array([-0., -0.]),
     'nJnet': array([-0., -0.]),
     'wFlux': array([2.38525e+14, 5.17849e+13]),
     'sFlux': array([2.38585e+14, 5.18003e+13]),
     'eFlux': array([2.38474e+14, 5.17765e+13]),
     'nFlux': array([2.38465e+14, 5.17969e+13]),
     'wDF': array([1.00029 , 0.970428]),
     'sDF': array([1.00054 , 0.970716]),
     'eDF': array([1.00008 , 0.970271]),
     'nDF': array([1.00004 , 0.970653]),
     'nwDF': array([1.00286 , 0.940338]),
     'neDF': array([1.00338 , 0.939482]),
     'seDF': array([1.00237 , 0.939798]),
     'swDF': array([1.00288 , 0.937232])}



Part 2 - get form factors
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    # Get serpent power density
    npins = 17
    vol_pin = 21.42**2 / npins**2
    detFile_455 = './serpent/Fuel_only/Fuel_455_2g_900K_det0.m'
    det_455 = serpentTools.read(detFile_455)
    fastPower_455 = det_455.detectors['power_fast'].tallies[0:npins, 0:npins]
    thermalPower_455 = det_455.detectors['power_thermal'].tallies[0:npins, 0:npins]
    serpentPower_455 = fastPower_455 + thermalPower_455
    xs_455 = AFEN.getFormFactors(xs=xs_455, power_fast=fastPower_455, power_thermal=thermalPower_455, vol=vol_pin)

    detFile_260 = './serpent/Fuel_only/Fuel_260_2g_det0.m'
    det_260 = serpentTools.read(detFile_260)
    fastPower_260 = det_260.detectors['power_fast'].tallies[0:npins, 0:npins]
    thermalPower_260 = det_260.detectors['power_thermal'].tallies[0:npins, 0:npins]
    serpentPower_260 = fastPower_260 + thermalPower_260
    xs_260 = AFEN.getFormFactors(xs=xs_260, power_fast=fastPower_260, power_thermal=thermalPower_260, vol=vol_pin)

    # fig, axs = plt.subplots(2,2)

    # form factors for 4.55% assembly
    AFEN.meshPlot(res2d=xs_455['ff'][0], npins=17, cbar_label='FF', xlabel='X (cm)', ylabel='Y (cm)', fontsize=13)
    AFEN.meshPlot(res2d=xs_455['ff'][1], npins=17, cbar_label='FF', xlabel='X (cm)', ylabel='Y (cm)', fontsize=13)

    # form factors for 2.60% assembly
    AFEN.meshPlot(res2d=xs_260['ff'][0], npins=17, cbar_label='FF', xlabel='X (cm)', ylabel='Y (cm)', fontsize=13)
    AFEN.meshPlot(res2d=xs_260['ff'][1], npins=17, cbar_label='FF', xlabel='X (cm)', ylabel='Y (cm)', fontsize=13)


.. parsed-literal::

    kappaSigmaF*Flx*Vol (homogenous) = 45881.84688206413
    kappaSigmaF*Flx*Vol (heterogenous) = 45881.640158760005
    kappaSigmaF*Flx*Vol (homogenous) = 45881.79101492478
    kappaSigmaF*Flx*Vol (heterogenous) = 45881.64174636001



.. image:: hw6_files/final_project_12_1.png



.. image:: hw6_files/final_project_12_2.png



.. image:: hw6_files/final_project_12_3.png



.. image:: hw6_files/final_project_12_4.png


Part 3 - Use the NEM solver to get the surface homogenous fluxes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    # Now set the entire problem up. Solve NEM in x and y directions
    resFile = './serpent/SMR/SMR_Ref_2D_2g_res.m'
    _, bc_f12 =  NEM.GetSerpentRes(resFile, 'F12', 0)
    xs_ref, bc_ref =  NEM.GetSerpentRes(resFile, 'Ref', 0)
    _, bc_f2 =  NEM.GetSerpentRes(resFile, 'F2', 0)
    _, bc_f11 =  NEM.GetSerpentRes(resFile, 'F11', 0)

    #########################################
    # Setup transverse leakages for each node
    #########################################
    dx, dy = 21.42, 21.42  # cm  - assembly length

    # F12
    trLeakage_f12 = {}
    trLeakage_f12['eL'] = bc_ref['nJnet'] - bc_ref['sJnet']
    trLeakage_f12['eD'] = xs_ref['diff'] # use normal diff here.
    trLeakage_f12['edx'] = dx

    trLeakage_f12['wL'] = bc_f12['nJnet'] - bc_f12['sJnet'] # use values from the node itself (since it is reflective ...)
    trLeakage_f12['wD'] = xs_260['cmmdiff']
    trLeakage_f12['wdx'] = dx

    trLeakage_f12['nL'] = bc_f12['eJnet'] - bc_f12['wJnet']
    trLeakage_f12['nD'] = xs_260['cmmdiff']
    trLeakage_f12['ndy'] = dy

    trLeakage_f12['sL'] = bc_f2['eJnet'] - bc_f2['wJnet']
    trLeakage_f12['sD'] = xs_455['cmmdiff']
    trLeakage_f12['sdy'] = dy

    # F2
    trLeakage_f2 = {}
    trLeakage_f2['eL'] = bc_f11['nJnet'] - bc_f11['sJnet']
    trLeakage_f2['eD'] = xs_260['cmmdiff']
    trLeakage_f2['edx'] = dx

    trLeakage_f2['wL'] = bc_f2['nJnet'] - bc_f2['sJnet'] # use values from the node itself (since it is reflective ...)
    trLeakage_f2['wD'] = xs_455['cmmdiff']
    trLeakage_f2['wdx'] = dx

    trLeakage_f2['nL'] = bc_f12['eJnet'] - bc_f12['wJnet']
    trLeakage_f2['nD'] = xs_260['cmmdiff']
    trLeakage_f2['ndy'] = dy

    trLeakage_f2['sL'] = bc_f2['eJnet'] - bc_f2['wJnet']
    trLeakage_f2['sD'] = xs_455['cmmdiff']
    trLeakage_f2['sdy'] = dy

    # Ref
    trLeakage_ref = {}
    trLeakage_ref['eL'] = bc_ref['nJnet'] - bc_ref['sJnet']
    trLeakage_ref['eD'] = xs_ref['diff']
    trLeakage_ref['edx'] = dx

    trLeakage_ref['wL'] = bc_f12['nJnet'] - bc_f12['sJnet'] # use values from the node itself (since it is reflective ...)
    trLeakage_ref['wD'] = xs_260['cmmdiff']
    trLeakage_ref['wdx'] = dx

    trLeakage_ref['nL'] = bc_ref['eJnet'] - bc_ref['wJnet']
    trLeakage_ref['nD'] = xs_ref['diff']
    trLeakage_ref['ndy'] = dy

    trLeakage_ref['sL'] = bc_f11['eJnet'] - bc_f11['wJnet']
    trLeakage_ref['sD'] = xs_260['cmmdiff']
    trLeakage_ref['sdy'] = dy

    # F11
    trLeakage_f11 = {}
    trLeakage_f11['eL'] = bc_f11['nJnet'] - bc_f11['sJnet']
    trLeakage_f11['eD'] = xs_260['cmmdiff']
    trLeakage_f11['edx'] = dx

    trLeakage_f11['wL'] = bc_f2['nJnet'] - bc_f2['sJnet'] # use values from the node itself (since it is reflective ...)
    trLeakage_f11['wD'] = xs_455['cmmdiff']
    trLeakage_f11['wdx'] = dx

    trLeakage_f11['nL'] = bc_ref['eJnet'] - bc_ref['wJnet']
    trLeakage_f11['nD'] = xs_ref['diff']
    trLeakage_f11['ndy'] = dy

    trLeakage_f11['sL'] = bc_f11['eJnet'] - bc_f11['wJnet']
    trLeakage_f11['sD'] = xs_260['cmmdiff']
    trLeakage_f11['sdy'] = dy



    # Setup bcs for next step (get from NEM output)
    nem_flux = {}
    nem_flux['F12'] = {}
    nem_flux['F2'] = {}
    nem_flux['ref'] = {}
    nem_flux['F11'] = {}

    # Solve for each in dx direction
    start_time = time.time()

    nem_f12 = NEM.CartesianNem1D(dx, dy, xs_260, bc_f12, trLeakage_f12, 'x')
    nem_f12.TransverseLeakageCoef()  # obtain the coefficients of the TL
    nem_f12.GetExpansionCoeffs()

    nem_f2 = NEM.CartesianNem1D(dx, dy, xs_455, bc_f2, trLeakage_f2, 'x')
    nem_f2.TransverseLeakageCoef()  # obtain the coefficients of the TL
    nem_f2.GetExpansionCoeffs()

    nem_ref = NEM.CartesianNem1D(dx, dy, xs_ref, bc_ref, trLeakage_ref, 'x')
    nem_ref.TransverseLeakageCoef()  # obtain the coefficients of the TL
    nem_ref.GetExpansionCoeffs()

    nem_f11 = NEM.CartesianNem1D(dx, dy, xs_260, bc_f11, trLeakage_f11, 'x')
    nem_f11.TransverseLeakageCoef()  # obtain the coefficients of the TL
    nem_f11.GetExpansionCoeffs()
    end_time = time.time()

    nem_flux['F12']['w'] = nem_f12.westHomFlux
    nem_flux['F2']['w'] = nem_f2.westHomFlux
    nem_flux['ref']['w'] = nem_ref.westHomFlux
    nem_flux['F11']['w'] = nem_f11.westHomFlux

    nem_flux['F12']['e'] = nem_f12.eastHomFlux
    nem_flux['F2']['e'] = nem_f2.eastHomFlux
    nem_flux['ref']['e'] = nem_ref.eastHomFlux
    nem_flux['F11']['e'] = nem_f11.eastHomFlux

    print("NEM x direction is done, time is", end_time-start_time)

    # Solve for y
    start_time = time.time()
    nem_f12 = NEM.CartesianNem1D(dx, dy, xs_260, bc_f12, trLeakage_f12, 'y')
    nem_f12.TransverseLeakageCoef()  # obtain the coefficients of the TL
    nem_f12.GetExpansionCoeffs()

    nem_f2 = NEM.CartesianNem1D(dx, dy, xs_455, bc_f2, trLeakage_f2, 'y')
    nem_f2.TransverseLeakageCoef()  # obtain the coefficients of the TL
    nem_f2.GetExpansionCoeffs()

    nem_ref = NEM.CartesianNem1D(dx, dy, xs_ref, bc_ref, trLeakage_ref, 'y')
    nem_ref.TransverseLeakageCoef()  # obtain the coefficients of the TL
    nem_ref.GetExpansionCoeffs()

    nem_f11 = NEM.CartesianNem1D(dx, dy, xs_260, bc_f11, trLeakage_f11, 'y')
    nem_f11.TransverseLeakageCoef()  # obtain the coefficients of the TL
    nem_f11.GetExpansionCoeffs()
    end_time = time.time()
    print("NEM y direction is done, time is", end_time-start_time)

    nem_flux['F12']['n'] = nem_f12.northHomFlux
    nem_flux['F2']['n'] = nem_f2.northHomFlux
    nem_flux['ref']['n'] = nem_ref.northHomFlux
    nem_flux['F11']['n'] = nem_f11.northHomFlux

    nem_flux['F12']['s'] = nem_f12.southHomFlux
    nem_flux['F2']['s'] = nem_f2.southHomFlux
    nem_flux['ref']['s'] = nem_ref.southHomFlux
    nem_flux['F11']['s'] = nem_f11.southHomFlux

    for key in nem_flux.keys():
      for key2 in nem_flux[key].keys():
        print(key, "|", key2, "=", nem_flux[key][key2])


.. parsed-literal::

    NEM x direction is done, time is 0.14284634590148926
    NEM y direction is done, time is 0.03863716125488281
    F12 | w = [2.52434290e+14 5.69928506e+13]
    F12 | e = [1.80608491e+14 2.10505844e+13]
    F12 | n = [2.15397734e+14 4.72887709e+13]
    F12 | s = [2.53214832e+14 4.52170421e+13]
    F2 | w = [2.87056623e+14 3.80176960e+13]
    F2 | e = [2.48830002e+14 4.19079191e+13]
    F2 | n = [2.48827575e+14 4.19075340e+13]
    F2 | s = [2.87058751e+14 3.80180183e+13]
    ref | w = [1.87221735e+14 2.51751224e+13]
    ref | e = [1.17016584e+14 1.31356617e+12]
    ref | n = [1.17017310e+14 1.31379707e+12]
    ref | s = [1.87220887e+14 2.51747557e+13]
    F11 | w = [2.53189608e+14 4.51887577e+13]
    F11 | e = [2.15413081e+14 4.72931002e+13]
    F11 | n = [1.80645325e+14 2.10596265e+13]
    F11 | s = [2.52414891e+14 5.69849894e+13]


Part 4 - Use het. fluxes from the reference colorset calc + corner DF’s to obtain homogenous corner fluxes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    # DF's stored in bc_455 (F2) and bc_260 (F12)
    # corner fluxes stored in bc_fNNN (bc_f11,f12,f2,ref)
    # do 1/N * sum(flx_het_k

    # assemble serpent fluxes



    # corner point fluxes
    def getCornerFlux(asm):
      """
      Solves equation 13.2 for each corner.
      asm = Serpent heterogenous fluxes
      """
      flux_c = {}

      # Equation 13.2
      flux_c['nw'] = -asm['flux'] + asm['nFlux'] + asm['wFlux']
      flux_c['ne'] = -asm['flux'] + asm['nFlux'] + asm['eFlux']
      flux_c['se'] = -asm['flux'] + asm['sFlux'] + asm['eFlux']
      flux_c['sw'] = -asm['flux'] + asm['sFlux'] + asm['wFlux']
      flux_c['w'] = asm['wFlux']
      flux_c['e'] = asm['eFlux']
      flux_c['s'] = asm['sFlux']
      flux_c['n'] = asm['nFlux']
      flux_c['av'] = asm['flux']

      return flux_c

    # Using Equation 13.2
    flux_c_f12 = getCornerFlux(asm=bc_f12)
    flux_c_ref = getCornerFlux(asm=bc_ref)
    flux_c_f2 = getCornerFlux(asm=bc_f2)
    flux_c_f11 = getCornerFlux(asm=bc_f11)

    # Using Equation 13.3
    flux_bar_c_f12 = flux_c_f12
    flux_bar_c_f12['nw'] = flux_c_f12['nw']
    flux_bar_c_f12['ne'] = flux_c_f12['ne']
    flux_bar_c_f12['se'] = (flux_c_f12['se']*bc_260['seDF'] + flux_c_f2['ne']*bc_455['neDF'] + flux_c_f11['nw']*bc_260['nwDF'])/3/bc_260['seDF']
    flux_bar_c_f12['sw'] = (flux_c_f12['sw']*bc_260['swDF'] + flux_c_f2['nw']*bc_455['nwDF'] ) / 2 / bc_260['swDF']

    flux_bar_c_f2 = flux_c_f2
    flux_bar_c_f2['nw'] = (flux_c_f2['nw']*bc_455['nwDF'] + flux_c_f12['sw']*bc_260['swDF'])/2/bc_455['nwDF']
    flux_bar_c_f2['ne'] = (flux_c_f2['ne']*bc_455['neDF'] + flux_c_f12['se']*bc_260['seDF'] + flux_c_f11['nw']*bc_260['nwDF'])/3/bc_455['neDF']
    flux_bar_c_f2['se'] = (flux_c_f2['se']*bc_455['seDF'] + flux_c_f11['sw']*bc_260['swDF'])/2/bc_455['seDF']
    flux_bar_c_f2['sw'] = flux_c_f2['sw']

    flux_bar_c_f11 = flux_c_f11
    flux_bar_c_f11['nw'] = (flux_c_f11['nw']*bc_260['nwDF'] + flux_c_f2['ne']*bc_455['neDF'] + flux_c_f12['se']*bc_260['seDF'])/3/bc_260['nwDF']
    flux_bar_c_f11['ne'] = flux_c_f11['ne']
    flux_bar_c_f11['se'] = flux_c_f11['se']
    flux_bar_c_f11['sw'] = (flux_c_f11['sw']*bc_260['swDF'] + flux_c_f2['se']*bc_455['seDF'])/2/bc_260['swDF']


.. code:: ipython3


    g = 0
    flx = flux_bar_c_f12
    print("01: ", flx['nw'])
    print("11: ", flx['ne'])
    print("10: ", flx['se'])
    print("00: ", flx['sw'])
    print("y1: ", flx['n'])
    print("x1: ", flx['e'])
    print("y0: ", flx['s'])
    print("x0: ", flx['w'])
    print("phi_bar", flx['av'])


.. parsed-literal::

    01:  [2.42907e+14 5.04895e+13]
    11:  [1.71540e+14 1.93302e+13]
    10:  [2.17414547e+14 2.88977023e+13]
    00:  [2.77655747e+14 4.93364437e+13]
    y1:  [2.16384e+14 4.59049e+13]
    x1:  [1.83113e+14 2.21544e+13]
    y0:  [2.52849e+14 4.32762e+13]
    x0:  [2.54480e+14 5.33137e+13]
    phi_bar [2.27957e+14 4.87291e+13]


Question 4
~~~~~~~~~~

.. code:: ipython3

    # First get serpent power
    detFile = './serpent/SMR/SMR_Ref_2D_2g_det0.m'
    det = serpentTools.read(detFile)

    fastHetFlux_f2 = det.detectors['flux_fast'].tallies[0:npins, 0:npins]
    thermalHetFlux_f2 = det.detectors['flux_thermal'].tallies[0:npins, 0:npins]
    fastPower_f2 = det.detectors['power_fast'].tallies[0:npins, 0:npins]
    thermalPower_f2 = det.detectors['power_thermal'].tallies[0:npins, 0:npins]
    serpentPower_f2 = fastPower_f2 + thermalPower_f2

    fastHetFlux_f11 = det.detectors['flux_fast'].tallies[0:npins, npins:]
    thermalHetFlux_f11 = det.detectors['flux_thermal'].tallies[0:npins, npins:]
    fastPower_f11 = det.detectors['power_fast'].tallies[0:npins, npins:]
    thermalPower_f11 = det.detectors['power_thermal'].tallies[0:npins, npins:]
    serpentPower_f11 = fastPower_f11 + thermalPower_f11

    fastHetFlux_f12 = det.detectors['flux_fast'].tallies[npins:, 0:npins]
    thermalHetFlux_f12 = det.detectors['flux_thermal'].tallies[npins:, 0:npins]
    fastPower_f12 = det.detectors['power_fast'].tallies[npins:, 0:npins]
    thermalPower_f12 = det.detectors['power_thermal'].tallies[npins:, 0:npins]
    serpentPower_f12 = fastPower_f12 + thermalPower_f12

    fastHetFlux_ref = det.detectors['flux_fast'].tallies[npins:, npins:]
    thermalHetFlux_ref = det.detectors['flux_thermal'].tallies[npins:, npins:]
    fastPower_ref = det.detectors['power_fast'].tallies[npins:, npins:]
    thermalPower_ref = det.detectors['power_thermal'].tallies[npins:, npins:]

.. code:: ipython3

    # Now solving AFEM
    yvals = np.linspace(-dx/2, +dx/2, npins+1)
    yvals = 0.5*(yvals[1:]+yvals[0:-1])

    afen_f2 = AFEN.AFEN2D(xs_455, flux_bar_c_f2, dx, symbolic=True)
    afen_f2.ReconstructFlux()  # a built-in method to obtain the coeffs
    afen_f2.GetFlux2D(yvals, yvals)
    power_f2 = afen_f2.PinPower(ff=np.array([xs_455['ff'][0], xs_455['ff'][1]]))

    afen_f12 = AFEN.AFEN2D(xs_260, flux_bar_c_f12, dx, symbolic=True)
    afen_f12.ReconstructFlux()  # a built-in method to obtain the coeffs
    afen_f12.GetFlux2D(yvals, yvals)
    power_f12 = afen_f12.PinPower(ff=np.array([xs_260['ff'][0], xs_260['ff'][1]]))

    afen_f11 = AFEN.AFEN2D(xs_260, flux_bar_c_f11, dx, symbolic=True)
    afen_f11.ReconstructFlux()  # a built-in method to obtain the coeffs
    afen_f11.GetFlux2D(yvals, yvals)
    power_f11 = afen_f11.PinPower(ff=np.array([xs_260['ff'][0], xs_260['ff'][1]]))


    # ERROR
    error_f12 = (serpentPower_f12-power_f12)/serpentPower_f12*100
    error_f2 = (serpentPower_f2-power_f2)/serpentPower_f2*100
    error_f11 = (serpentPower_f11-power_f11)/serpentPower_f11*100

    fig, axs = plt.subplots(2,2, dpi=600)
    AFEN.meshPlot(error_f12, 17, fontsize=6, ax=axs[0][0], cbar_limits=(-10, 10), do_cbar=False)
    AFEN.meshPlot(error_f2, 17, fontsize=6, ax=axs[1][0], cbar_limits=(-10, 10), do_cbar=False)
    AFEN.meshPlot(error_f11, 17, fontsize=6, ax=axs[1][1], cbar_limits=(-10, 10), cbar_label='Percent Diff. Power')
    axs[0,1].remove()

    # AFEN
    fig, axs = plt.subplots(2,2, dpi=600)
    AFEN.meshPlot(power_f12, 17, fontsize=6, ax=axs[0][0], cbar_limits=(25, 175), do_cbar=False, cmap='inferno')
    AFEN.meshPlot(power_f2, 17, fontsize=6, ax=axs[1][0], cbar_limits=(25, 175), do_cbar=False, cmap='inferno')
    AFEN.meshPlot(power_f11, 17, fontsize=6, ax=axs[1][1], cbar_limits=(25, 175), cbar_label='Power (AFEN)', cmap='inferno')
    axs[0,1].remove()

    # REFERENCE
    fig, axs = plt.subplots(2,2, dpi=600)
    AFEN.meshPlot(serpentPower_f12, 17, fontsize=6, ax=axs[0][0], cbar_limits=(25, 175), do_cbar=False, cmap='inferno')
    AFEN.meshPlot(serpentPower_f2, 17, fontsize=6, ax=axs[1][0], cbar_limits=(25, 175), do_cbar=False, cmap='inferno')
    AFEN.meshPlot(serpentPower_f11, 17, fontsize=6, ax=axs[1][1], cbar_limits=(25, 175), cbar_label='Power (Serpent)', cmap='inferno')
    axs[0,1].remove()


.. parsed-literal::

    /tmp/ipykernel_409546/2576880132.py:22: RuntimeWarning: invalid value encountered in divide
      error_f12 = (serpentPower_f12-power_f12)/serpentPower_f12*100
    /tmp/ipykernel_409546/2576880132.py:23: RuntimeWarning: invalid value encountered in divide
      error_f2 = (serpentPower_f2-power_f2)/serpentPower_f2*100
    /tmp/ipykernel_409546/2576880132.py:24: RuntimeWarning: invalid value encountered in divide
      error_f11 = (serpentPower_f11-power_f11)/serpentPower_f11*100



.. image:: hw6_files/final_project_20_1.png



.. image:: hw6_files/final_project_20_2.png



.. image:: hw6_files/final_project_20_3.png


.. code:: ipython3

    def get_mean_abs(arr):
      summed = 0.0
      count = int(0)
      this_max = 0.0
      for row in arr:
        for col in row:
          if np.isnan(col):
            pass
          else:
            count += 1
            summed += np.abs(col)
            if np.abs(col) > this_max:
              this_max = np.abs(col)
      return summed/count, this_max

    print("mean abs f12", get_mean_abs(error_f12))
    print("mean abs f11", get_mean_abs(error_f11))
    print("mean abs f2", get_mean_abs(error_f2))


.. parsed-literal::

    mean abs f12 (np.float64(1.8525427853684189), np.float64(7.740551568722151))
    mean abs f11 (np.float64(1.7959443257000267), np.float64(5.945815885279578))
    mean abs f2 (np.float64(1.975275112565977), np.float64(8.601424492093734))

