# -*- coding: utf-8 -*-
"""
Module:

pointsource_sphere

Description of this file:
    - Solution of the point source within a sphere
    - The heterogenous sphere is divided into multiple shells
    - Only the absorption reaction is considered here for simplicity
    - It is assumed that the source does not absrob the neutrons
    - Although the source radius is defined, it is over-written with zero

Created on Sat May 17 18:20:00 2025 @author: Dan Kotlyar
Last updated on Mon May 19 11:30:00 2025 @author: Dan Kotlyar

email: dan.kotlyar@me.gatech.edu


"""

import time
import matplotlib.pyplot as plt
# from matplotlib import rcParams
import numpy as np
# rcParams['figure.dpi'] = 300

FONT_SIZE = 16
MARKER_SIZE = 6

ALLOWED_SCHEMES = ['st', 'dt', 'analytic']  # surface / delta tracking


class WeightedNeutron:
    """
    A class for neutrons. Kinda simple right now but we
    would ideally use this more and more and expand on it as
    more methods are built/added.

    Only currently used in the weighted delta tracking routine.

    Parameters
    ----------

    Attributes
    ----------
    inps : weight
        weight of the neutron

    """
    def __init__(self):
        """
        """
        self.weight = None # Set to None to ensure that we get an error if we forget to set starting weight

    def setWeight(self, weight: float):
        """
        Updates the weight of the neutron.

        Parameters
        ----------
        weight : float
            new weight of the neutron

        """
        if (weight >= 0.0) & (weight  <= 1.0):
            self.weight = weight
        else:
            raise ValueError("Weight value must be between 0 and 1.0")


class PointSourceInSphere():
    """MC solution for a point source positioned in an heterogenous sphere

    Parameters
    ----------
    nMC : int
        number of repetitive monte carlo simulations
    S0 : int
        number of source neutrons
    R : float
        radius of the sphere in cm
    sigT : array
        total cross section array for all the different regions


    Attributes
    ----------
    inps : dict
        container with inputs
    resAN : dict
        container with analytic results
    resST : dict
        container with surface tracking results
    resDT : dict
        container with delta tracking results
    times : dict
        container with execution times

    Raises
    -------
    TypeError
        if nMC is not int
    ValueError
        if nMC is not positive

    """

    def __init__(self, nMC: int, S0: int, R: float, sigT: np.ndarray | list):

        # check errors (example)
        if not isinstance(nMC, int):
            raise TypeError("nMC must be integer and not {}".format(nMC))
        if nMC <= 0:
            raise ValueError("nMC must be positive int and not {}".format(nMC))

        self.resST = {}  # surface tracking results
        self.resDT = {}  # delta tracking
        self.resWDT = {} # weighted delta tracking results
        self.resAN = {}  # analytic results
        self.times = {}  # solution times

        nreg = len(sigT)
        vol = (4/3) * np.pi * R**3 / nreg  # single volume shell
        # define equal-volume radii (all shells are divided into equal vol which means different R)
        radii = np.zeros(nreg)
        radii[0] = (vol*3/(4*np.pi))**(1/3)
        for i in range(1, nreg):
            radii[i] = (radii[i-1]**3 + vol*3/(4*np.pi))**(1/3)

        # save inputs dict
        self.inps = {'R': R, 'nMC': nMC, 'sigT': sigT, 'radii': radii,
                     'nreg': nreg, 'S0': S0}



    def Solve(self, scheme):
        """
        MC solution using a specific random walk scheme

        Parameters
        ----------
        scheme : str
            neutron tracking scheme ['st', 'dt', 'analytic', 'wdt']

        """

        nMC = self.inps['nMC']
        nreg = self.inps['nreg']

        if scheme == 'analytic':
            print("\n\n\nNow solving with method:", scheme)
            self._Analytic()

        elif scheme.lower() == 'st':
            print("Now solving with method:", scheme)
            start_time = time.time()
            flux = np.zeros((nreg, nMC))
            leakage = np.zeros(nMC)
            for i in range(nMC):
                if i % 10 == 0:
                    print("Running generation", i)
                flux[:, i], leakage[i] = self._SolveST()
            self.resST = {'flx': np.mean(flux, 1),
                          'errflx': np.std(flux, 1),
                          'leakage': np.mean(leakage),
                          'errleakage': np.std(leakage)}
            self.times['ST'] = time.time() - start_time
            print("Finished after", i, "MC runs! Time was", self.times['ST'])
            print("\n\n")

        elif scheme.lower() == 'dt':
            print("Now solving with method:", scheme)
            start_time = time.time()
            flux = np.zeros((nreg, nMC))
            leakage = np.zeros(nMC)
            for i in range(nMC):
                if i % 10 == 0:
                    print("Running generation", i)
                flux[:, i], leakage[i] = self._SolveDT()
            self.resDT = {'flx': np.mean(flux, 1),
                          'errflx': np.std(flux, 1),
                          'leakage': np.mean(leakage),
                          'errleakage': np.std(leakage)}
            self.times['DT'] = time.time() - start_time
            print("Finished after", i, "MC runs! Time was", self.times['DT'])
            print("\n\n")

        elif scheme.lower() == 'wdt':
            print("Now solving with method:", scheme)
            start_time = time.time()
            flux = np.zeros((nreg, nMC))
            leakage = np.zeros(nMC)
            for i in range(nMC):
                if i % 10 == 0:
                    print("Running generation", i)
                flux[:, i], leakage[i] = self._SolveWDT()
            self.resWDT = {'flx': np.mean(flux, 1),
                          'errflx': np.std(flux, 1),
                          'leakage': np.mean(leakage),
                          'errleakage': np.std(leakage)}
            self.times['WDT'] = time.time() - start_time
            print("Finished after", i, "MC runs! Time was", self.times['WDT'])
            print("\n\n")

        else:
            raise ValueError('scheme {} does not exist in {}'
                             .format(scheme,ALLOWED_SCHEMES))


    def _Analytic(self):
        """Analytic solution"""

        nreg = self.inps['nreg']
        radii = self.inps['radii']
        sigT = self.inps['sigT']

        # loop over all regions
        I0 = self.inps['S0']
        Ix_I0 = np.full(nreg, I0, dtype=float) # ix and i0 are the intensity.
        for j in range(nreg):

            if j > 0:
                dr = radii[j] - radii[j-1]
                Ix_I0[j] = Ix_I0[j-1] * np.exp(-sigT[j] * dr)
            else: # j == 0
                dr = radii[j]
                Ix_I0[j] = I0 * np.exp(-sigT[j] * dr)

        leakage = Ix_I0[-1] / I0 # last value over initial source = leakage
        self.resAN = {'flx': Ix_I0, 'leakage': leakage}


    def _SolveST(self):
        """MC solution using ray tracking"""

        nreg = self.inps['nreg']
        radii = self.inps['radii']
        sigT = self.inps['sigT']
        S0 = self.inps['S0']

        # reset neutron flux
        flux = np.zeros(nreg) # surface fluxes

        # sample particles position (R, theta, phi)
        xipos = np.random.rand(S0, 3) # random number psi

        # Particles are "born" in the most-inner sphere (R,teta,phi)
        # NOTE we are actualyl sampling a sphere but we hardcode it to be at r=0.0
        Rin = 0.0  # fake radius (just to show how a position is sampled)
        R0 = Rin * xipos[:,0]**(1/3)        # sample radius
        theta = np.arccos(1-2*xipos[:,1]) # teta angle
        phi = 2*np.pi*xipos[:,2]           # phi angle

        # Sample directional flight of neutrons
        xidir = np.random.rand(S0, 2) # xidir !!! not xipos
        wtheta = np.arccos(1-2*xidir[:,0])  # teta
        wphi = 2*np.pi*xidir[:,1]            # phi angle

        for i in range(S0): # Iterate through particles.

            # (X0,Y0,Z0) - before collision
            x0 = R0[i] * np.sin(theta[i]) * np.cos(phi[i])
            y0 = R0[i] * np.sin(theta[i]) * np.sin(phi[i])
            z0 = R0[i] * np.cos(theta[i])

            # (WX,WY,WZ) - direction (use wtheta and wphi instead of theta and phi)
            wx = np.sin(wtheta[i]) * np.cos(wphi[i])
            wy = np.sin(wtheta[i]) * np.sin(wphi[i])
            wz = np.cos(wtheta[i])

            # Solve the following inequety and find the intersection surfaces
            # Lx^2+Ly^2+Lz^2 = Ri^2 (eq.1)
            # Lx=x0+Si*wx ;Ly=y0+Si*wy; Lz=z0+Si*wz (eq.2)

            for j in range(nreg): # loop over all possible intersections
                # Sample collision's position
                xi = np.random.rand() # random number from Uniform [0,1]
                Si = -np.log(xi) / sigT[j] # Si is what we call L in the lecture notes.

                x1=x0 + Si*wx
                y1=y0 + Si*wy
                z1=z0 + Si*wz

                if  (x1**2 + y1**2 + z1**2) < radii[j]**2:
                    # neutron is absorbed in this shell
                    break
                elif j == nreg-1:
                    # particle escaped from the system
                    dist = radii[j] - radii[j-1]
                    flux[j] += 1  # score the neutron
                    break

                else:
                    # Dist is the distance travelled ebfore crossing a surface
                    if j > 0: # NOT First shell
                        dist = radii[j] - radii[j-1] # dist is diff between radii
                    else: # Just the first shell
                        dist = radii[j] # just distance from 0.0 to radii of first shell
                    flux[j] += 1  # score the neutron
                    # update the position of the neutron
                    x0=x0+dist*wx
                    y0=y0+dist*wy
                    z0=z0+dist*wz

        # get relative leakage fraction
        leakage = flux[-1] / S0

        return flux, leakage

    def _SolveDT(self):
        """MC solution using delta tracking"""

        nreg = self.inps['nreg']
        radii = self.inps['radii']
        sigT = self.inps['sigT']
        S0 = self.inps['S0']
        radii2 = radii**2  # radii squared
        # majorant cross section
        sigMaj = max(sigT)

        # reset neutron flux
        flux = np.zeros(nreg)

        # sample particles position (R, theta, phi)
        xipos = np.random.rand(S0, 3)

        # Particles are "born" in the most-inner sphere (R,teta,phi)
        Rin = 0.0  # fake radius (just to show how a position is sampled)
        R0 = Rin * xipos[:,0]**(1/3)        # sample radius
        theta = np.arccos( 1-2*xipos[:, 1] )  # teta angle
        phi = 2*np.pi*xipos[:, 2]             # phi angle

        # Sample directional flight of neutrons
        xidir = np.random.rand(S0, 2) # xidir !!! not xipos
        wtheta = np.arccos(1-2*xidir[:,0])  # teta
        wphi = 2*np.pi*xidir[:,1]            # phi angle

        for i in range(S0):

            # (X0,Y0,Z0) - before collision
            x0 = R0[i] * np.sin(theta[i]) * np.cos(phi[i])
            y0 = R0[i] * np.sin(theta[i]) * np.sin(phi[i])
            z0 = R0[i] * np.cos(theta[i])

            # (WX,WY,WZ) - direction
            wx = np.sin(wtheta[i]) * np.cos(wphi[i])
            wy = np.sin(wtheta[i]) * np.sin(wphi[i])
            wz = np.cos(wtheta[i])

            # this basically tells us what shell the neutron starts from ...
            idx1 = 0

            while 1:

                xi = np.random.rand() # random number from Uniform [0,1] #  sample xi so that we can sample a distance
                Si = -1.0 * np.log(xi) / sigMaj ####

                x1=x0 + Si*wx
                y1=y0 + Si*wy
                z1=z0 + Si*wz

                R2 = x1**2+y1**2+z1**2

                if R2 >= radii2[-1]:
                    # neutron leaked
                    # if it leaks it passes all the surfaces from idx1=0 up until the end
                    flux[:] += 1.0  # neutron crossed all shells so we score all the shells.
                    break
                else:
                    # neutron did not leak! -- did not leak so find new shell that it is in and update idx1
                    idx1 = np.where(R2 < radii2)[0][0] # returns index first case of True is found!
                    # woodcock rejection method
                    if np.random.rand() < sigT[idx1] / sigMaj:
                        # real collision
                        flux[0:idx1] += 1.0 # score up to but NOT INCLUDING INDEX 1
                        break
                    else:
                        # virtual collision
                        x0 = x1
                        y0 = y1
                        z0 = z1


        # get relative leakage fraction
        leakage = flux[-1] / S0

        return flux, leakage

    def _SolveWDT(self):
        """
        MC solution using weighted delta tracking routine.

        Uses a Russian Roulette threshold of 0.25 to discard particles.

        Call via PointSourceInSphere() object.

        Based on:

        L.W.G. Morgan, D. Kotlyar, 2015. "Weighted-delta-tracking for Monte Carlo particle transport,"
        Annals of Nuclear Energy, 85, 1184-1188.
        """

        # Get some information about the problem
        nreg = self.inps['nreg']
        radii = self.inps['radii']
        sigT = self.inps['sigT']
        S0 = self.inps['S0']
        radii2 = radii**2

        # Starting particle weight
        starting_weight = 1.0

        # Majorant XS
        sigMaj = max(sigT)

        # Reset neutron flux just in case!
        flux = np.zeros(nreg)

        # sample particles position (R, theta, phi)
        xipos = np.random.rand(S0, 3)

        # Particles are "born" in the most-inner sphere (R,teta,phi)
        Rin = 0.0  # fake radius (just to show how a position is sampled)
        R0 = Rin * xipos[:,0]**(1/3)        # sample radius
        theta = np.arccos( 1-2*xipos[:, 1] )  # teta angle
        phi = 2*np.pi*xipos[:, 2]             # phi angle

        # Sample directional flight of neutrons
        xidir = np.random.rand(S0, 2) # xidir !!! not xipos
        wtheta = np.arccos(1-2*xidir[:,0])  # teta
        wphi = 2*np.pi*xidir[:,1]            # phi angle

        # Now loop over all neutrons
        for i in range(S0):

            # Start with a beautiful neutron object and set its weight to starting weight
            this_neutron = WeightedNeutron()
            this_neutron.setWeight(weight=starting_weight)

            # (X0,Y0,Z0) - before collision
            x0 = R0[i] * np.sin(theta[i]) * np.cos(phi[i])
            y0 = R0[i] * np.sin(theta[i]) * np.sin(phi[i])
            z0 = R0[i] * np.cos(theta[i])

            # (WX,WY,WZ) - direction
            wx = np.sin(wtheta[i]) * np.cos(wphi[i])
            wy = np.sin(wtheta[i]) * np.sin(wphi[i])
            wz = np.cos(wtheta[i])

            # this basically tells us what shell the neutron starts from ...
            idx1 = 0

            # Loops until neutron is no longer with us :(
            while 1:

                xi = np.random.rand() # random number from Uniform [0,1] #  sample xi so that we can sample a distance
                Si = -1.0 * np.log(xi) / sigMaj ####

                # Collision distance
                x1=x0 + Si*wx
                y1=y0 + Si*wy
                z1=z0 + Si*wz

                R2 = x1**2+y1**2+z1**2

                if R2 >= radii2[-1]:
                    # neutron leaked
                    # if it leaks it passes all the surfaces from idx1=0 up until the end
                    flux[:] += this_neutron.weight  # neutron crossed all shells so we score all the shells.
                    break
                else:
                    # Neutron did not leak! -- did not leak so find new shell that it is in and update idx1
                    idx1 = np.where(R2 < radii2)[0][0] # returns index first case of True is found!

                    # woodcock rejection method
                    sigRatio_DT = sigT[idx1] / sigMaj # standard delta tracking eta value

                    # Get score that particle contributes and the particles new weight after the collision
                    score = sigRatio_DT * this_neutron.weight

                    # Update flux
                    flux[0:idx1] += score

                    # Update weight of neutron before continuing
                    new_weight = this_neutron.weight * (1.0 - sigRatio_DT)
                    this_neutron.setWeight(new_weight)

                    # Update initial location and continue.
                    x0 = x1
                    y0 = y1
                    z0 = z1

                    # Now play russian roullette to determine if we want to keep tracking
                    # Idea is to throw away low weight particles / regenerate daughter particles while ensuring a fair game
                    if this_neutron.weight < (0.25 / starting_weight):
                        # Roll for whether it survives.
                        p = np.random.rand()
                        if p < this_neutron.weight:
                            # We have survived russian roulette (phew!)
                            this_neutron.setWeight(weight=starting_weight)
                        else:
                            # particle is kill
                            break

        # get relative leakage fraction
        leakage = flux[-1] / S0

        return flux, leakage

    def PlotFluxes(self, scheme):
        """

        Plots neutron fluxes for a given method.

        Parameters
        ----------
        scheme : str
            scheme to plot fluxes for ('dt', 'wdt', 'all', 'st')

        """

        if not isinstance(scheme, str):
            raise TypeError("scheme must be str and not {}".format(scheme))

        radii = self.inps['radii']
        ylabel = "Flux [#]"
        xlabel = "Distance from center, cm"
        mfc = "white"  # marker fill "white" / None
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        if 'flx' in self.resAN:
            plt.plot(radii, self.resAN['flx'], 'k-', ms=MARKER_SIZE)
        if scheme.lower() == 'st':
            if 'flx' in self.resST:
                plt.errorbar(radii, self.resST['flx'],
                             yerr=self.resST['errflx'],
                             fmt='go', mfc=mfc, capsize=5)
        if scheme.lower() == 'dt':
            if 'flx' in self.resDT:
                plt.errorbar(radii, self.resDT['flx'],
                             yerr=self.resDT['errflx'],
                             fmt='r*', capsize=5)
        if scheme.lower() == 'wdt':
            if 'flx' in self.resWDT:
                plt.errorbar(radii, self.resWDT['flx'],
                             yerr=self.resWDT['errflx'],
                             fmt='bs', capsize=5,
                             markerfacecolor='w',
                             markersize=3)
        if scheme.lower() == 'all':
            # analytical
            plt.plot(radii, self.resAN['flx'], 'k-', ms=MARKER_SIZE,
                             label='Analytical')
            # surface tracking
            if 'flx' in self.resST:
                plt.errorbar(radii, self.resST['flx'],
                             yerr=self.resST['errflx'],
                             fmt='go', mfc=mfc, capsize=5,
                             label='ST')
            # delta tracking
            if 'flx' in self.resDT:
                plt.errorbar(radii, self.resDT['flx'],
                             yerr=self.resDT['errflx'],
                             fmt='r*', capsize=5,
                             label='DT')
            # weighted delta tracking
            if 'flx' in self.resWDT:
                plt.errorbar(radii, self.resWDT['flx'],
                             yerr=self.resWDT['errflx'],
                             fmt='bs', capsize=5,
                             markerfacecolor='w',
                             markersize=3,
                             label='WDT')

        plt.legend()

        plt.grid()
        plt.rc('font', size=FONT_SIZE)        # text sizes
        plt.rc('axes', labelsize=FONT_SIZE)   # labels
        plt.rc('xtick', labelsize=FONT_SIZE)  # tick labels
        plt.rc('ytick', labelsize=FONT_SIZE)  # tick labels

    def PlotDifferences(self):
        """
        Plots differences in fluxes
        """


        radii = self.inps['radii']
        ylabel = "Difference in Flux [%]"
        xlabel = "Distance from center, cm"
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)

        diffST = 100*(1-self.resST['flx']/self.resAN['flx'])
        diffDT = 100*(1-self.resDT['flx']/self.resAN['flx'])
        diffWDT = 100*(1-self.resWDT['flx']/self.resAN['flx'])

        plt.plot(radii, diffST, 'g--o', mfc='white', ms=MARKER_SIZE)
        plt.plot(radii, diffDT, 'r--*', ms=MARKER_SIZE)
        plt.plot(radii, diffWDT, 'b--s', ms=3, markerfacecolor='w')

        plt.legend(['ST', 'DT', 'WDT'])

        plt.grid()
        plt.rc('font', size=FONT_SIZE)        # text sizes
        plt.rc('axes', labelsize=FONT_SIZE)   # labels
        plt.rc('xtick', labelsize=FONT_SIZE)  # tick labels
        plt.rc('ytick', labelsize=FONT_SIZE)  # tick labels

    def PrintFOM(self, scheme: str):
        """
        Prints FOM for the desired scheme

        Uses the average of the variance squared and time for the FOM:
            FOM = 1/(sig^2 * time)
            where sig^2 = mean(sig^2(xyz)])

        Parameters
        ----------
        scheme : str
            scheme to plot fluxes for ('dt', 'wdt', 'all', 'st')


        """
        if scheme == 'analytic':
            print("No variance for analytical solution!")
        elif scheme.lower() == 'st':
            avg_dev = np.mean(self.resST['errflx']**2)
            time = self.times['ST']
        elif scheme.lower() == 'dt':
            avg_dev = np.mean(self.resDT['errflx']**2)
            time = self.times['DT']
        elif scheme.lower() == 'wdt':
            avg_dev = np.mean(self.resWDT['errflx']**2)
            time = self.times['WDT']
        else:
            raise ValueError('scheme {} does not exist in {}'
                             .format(scheme,ALLOWED_SCHEMES))
        FOM = 1.0 / avg_dev / time
        print("FOM - "+scheme+":", FOM)

def PlotSigT(mc: PointSourceInSphere, yLower: float, yUpper: float):
    """
    Def for plotting sigma_total.

    Parameters
    ----------
    mc :
        A PointSourceInSphere() object
    yLower : int
        Lower bound for y axis on plot
    yUpper : array
        Upper bound for y axis on plot

    """


    widths = mc.inps['radii'] - np.append([0], mc.inps['radii'])[0:-1]
    mids = ( mc.inps['radii'] + np.append([0], mc.inps['radii'])[0:-1] ) / 2
    vals = mc.inps['sigT']

    starting_x = 0.0
    starting_y = 0.0
    x = []
    y = []

    for idx, w in enumerate(widths):
        ending_x = starting_x + w
        ending_y = vals[idx]

        x += [starting_x, starting_x, ending_x]
        y += [starting_y, ending_y, ending_y]

        starting_x = ending_x
        starting_y = ending_y
    x += [x[-1]]
    y += [0.0]

    # plt.bar(mids, mc_c.inps['sigT'], width=widths*1.02, edgecolor='black', facecolor='white')
    plt.plot(x, y, 'k-', linewidth=2)
    plt.ylim([yLower, yUpper])
    plt.ylabel('$\sigma_t$ $(cm^{-1})$', fontsize=15)
    plt.xlabel('Radius (cm)', fontsize=15)
