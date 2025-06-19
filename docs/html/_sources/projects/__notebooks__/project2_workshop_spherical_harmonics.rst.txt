.. _proj2_sph_harm_notebook:

NRE7203: Advanced Reactor Physics
=================================

Copyright (c) Dan Kotlyar

Under Copyright law, you do not have the right to provide these notes to
anyone else or to make any commercial use of them without express prior
permission from me.

.. code:: ipython3

    from sph_harmonics import getSphHarm, plotSphHarm, animateSphHarm, BibicallyAccurateSphericalHarmonics, imageCropper, gifMaker
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import numpy as np
    from IPython.display import Image


.. code:: ipython3

    # Lets make a nice plot of our favorite spherical harmonic function. We use the 'Bibically Accurate' version to show the shape of the spherical harmonics rather than a boring sphere
    l=1
    m=1
    elev = 25
    azimRot = 45
    nangles = 100
    dpi=100

    # Plot and save
    BibicallyAccurateSphericalHarmonics(l=l, m=m, sphType='real', nangles=nangles, elev=elev, azimRot=azimRot, dpi=dpi, filename='L1_m1.png', doPlot=True)

    # Crop image for convenience.
    imageCropper(filename='L1_m1.png', outname='cropped.png', left_inches=0.5, right_inches=0.4, top_inches=1.2, bottom_inches=1.15, dpi=dpi)




.. parsed-literal::

    /home/jonathon/projects/arp_private/workshop_05_spherical_harmonics_vStudents/sph_harmonics.py:221: UserWarning: The figure layout has changed to tight
      plt.tight_layout()


.. parsed-literal::

    saving figure as L1_m1.png
    done saving figure ...



.. image:: project2_workshop_spherical_harmonics_files/project2_workshop_spherical_harmonics_3_2.png


.. code:: ipython3

    # Lets make a nice plot of our favorite spherical harmonic function. We use the 'Bibically Accurate' version to show the shape of the spherical harmonics rather than a boring sphere
    l=1
    m=0
    elev = 25
    azimRot = 45
    nangles = 100
    dpi=100

    # Plot and save
    BibicallyAccurateSphericalHarmonics(l=l, m=m, sphType='real', nangles=nangles, elev=elev, azimRot=azimRot, dpi=dpi, filename='L1_m0.png', doPlot=True)

    # Crop image for convenience.
    imageCropper(filename='L1_m0.png', outname='cropped.png', left_inches=0.5, right_inches=0.4, top_inches=1.2, bottom_inches=1.15, dpi=dpi)




.. parsed-literal::

    saving figure as L1_m0.png
    done saving figure ...



.. image:: project2_workshop_spherical_harmonics_files/project2_workshop_spherical_harmonics_4_1.png


.. code:: ipython3

    # Lets make a nice plot of our favorite spherical harmonic function. We use the 'Bibically Accurate' version to show the shape of the spherical harmonics rather than a boring sphere
    l=1
    m=-1
    elev = 25
    azimRot = 45
    nangles = 100
    dpi=100

    # Plot and save
    BibicallyAccurateSphericalHarmonics(l=l, m=m, sphType='real', nangles=nangles, elev=elev, azimRot=azimRot, dpi=dpi, filename='L1_m-1.png', doPlot=True)

    # Crop image for convenience.
    imageCropper(filename='L1_m-1.png', outname='cropped.png', left_inches=0.5, right_inches=0.4, top_inches=1.2, bottom_inches=1.15, dpi=dpi)




.. parsed-literal::

    saving figure as L1_m-1.png
    done saving figure ...



.. image:: project2_workshop_spherical_harmonics_files/project2_workshop_spherical_harmonics_5_1.png


Now lets make a cool gif!
=========================

Do not run this with large dpi or nFrames unless you value your computer’s memory and your time
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    # Lets make a gif (do not run unless you value your time)
    from IPython.display import Image

    l=1
    m=0
    nangles = 25
    dpi=50
    nFrames = 200

    elev_angles = np.linspace(0, 2*np.pi, nFrames) * 360 / 2.0 / np.pi # 1 full 360 degree rotation
    azim_angles = np.linspace(0, 2*np.pi, nFrames) * 360 / 2.0 / np.pi # 1 full 360 degree rotation


    filename_list = []
    for i, _ in enumerate(elev_angles):
      filename = './L1_M0_GIF/'+str(i)+'.png'
      BibicallyAccurateSphericalHarmonics(l=l, m=m, sphType='real', nangles=nangles, elev=elev_angles[i], azimRot=azim_angles[i], dpi=dpi, filename=filename, doPlot=False)
      imageCropper(filename=filename, outname=filename, left_inches=0.5, right_inches=0.4, top_inches=1.2, bottom_inches=1.15, dpi=dpi)
      filename_list.append(filename)
    gifMaker(output_filename='./L1_M0_GIF/L1_M0.gif', duration=5, filenames=filename_list)




.. parsed-literal::

    saving figure as ./L1_M0_GIF/0.png
    done saving figure ...
    saving figure as ./L1_M0_GIF/1.png
    done saving figure ...
    saving figure as ./L1_M0_GIF/2.png
    ...


.. code:: ipython3

    Image(url='./L1_M1_GIF/L1_M1.gif')




.. image:: project2_workshop_spherical_harmonics_files/L1_M1.gif
	:align: center
	:width: 2000



Inclass Stuff
-------------

.. code:: ipython3

    # Plot the real part
    plotSphHarm(l, m, sphType='real', elev=40)



.. image:: project2_workshop_spherical_harmonics_files/project2_workshop_spherical_harmonics_10_0.png


.. code:: ipython3

    # Plot the imaginary part
    plotSphHarm(l, m, sphType='imag', elev=40)


.. parsed-literal::

    /home/jonathon/projects/arp_private/workshop_05_spherical_harmonics_vStudents/sph_harmonics.py:68: RuntimeWarning: invalid value encountered in divide
      fcolors = (fcolors - fmin) / (fmax - fmin)



.. image:: project2_workshop_spherical_harmonics_files/project2_workshop_spherical_harmonics_11_1.png


.. code:: ipython3

    # animate
    nFrames = 5
    animateSphHarm(nFrames, l, m, sphType='real')



.. image:: project2_workshop_spherical_harmonics_files/project2_workshop_spherical_harmonics_12_0.png



.. image:: project2_workshop_spherical_harmonics_files/project2_workshop_spherical_harmonics_12_1.png



.. image:: project2_workshop_spherical_harmonics_files/project2_workshop_spherical_harmonics_12_2.png



.. image:: project2_workshop_spherical_harmonics_files/project2_workshop_spherical_harmonics_12_3.png



.. image:: project2_workshop_spherical_harmonics_files/project2_workshop_spherical_harmonics_12_4.png

