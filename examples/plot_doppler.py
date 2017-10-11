#!/usr/bin/env python

import sys
import matplotlib.pyplot as plt
import numpy as np

sys.path.append('..')
from keplerOrbit import KeplerOrbit


from orbparlib import readfile_OrbitalElements


if __name__=="__main__" :


    parfname = "orbelems/xper.param"
    par = readfile_OrbitalElements(parfname)
    ax    = par[8]
    porb  = par[4]
    pdot  =  par[5]
    ecc   =  par[11]
    epoch =  par[7]
    #Omega =  Omega
    periapse =  par[9]
    arate    =  par[10]
    #incl    =   incl, 
    epoch_type = par[6]

    
    plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    
    Omega=0.0
    incl =90.0
    
    mjd0_JD = 2400000.5 
    if mjd0_JD < epoch :
        epoch_mjd = epoch - mjd0_JD
    else :
        epoch_mjd = epoch 
    
    ko = KeplerOrbit()
    ko.setOrbitalElements2(ax, porb, pdot, ecc, epoch_mjd, Omega, periapse, arate, incl);
    if epoch_type==0 :
        ko.setEpochTypeP()
    else :
        ko.setEpochTypeT()
    
    mjd_start = 55050
    #mjd_end   = 57500
    mjd_end   = 58000
    mjd_step  = 1
    
    
    vmjd = np.arange(mjd_start, mjd_end, mjd_step)
    nv = len(vmjd)
    vdplr = np.zeros(nv)
    for idx in range(nv) :
        mjd = float(vmjd[idx])

        ###dplr = ko.getDopplerFactor(mjd)

        #dplr = ko.getDopplerFactor(mjd)-1
        dplr = (1./ko.getDopplerFactor(mjd) -1)*836  ### X Per period
        #dplr = (1./ko.getDopplerFactor(mjd) -1)

        vdplr[idx] = dplr
    
    plt.plot(vmjd, vdplr)
    #plt.axes().set_aspect('equal', 'datalim')
    
    plt.show()
