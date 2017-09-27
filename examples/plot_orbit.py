#!/usr/bin/env python

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

sys.path.append('..')
from keplerOrbit import KeplerOrbit

from orbparlib import readfile_OrbitalElements

#def get_orbitxyz(ax, porb, ecc, epoch, Omega, periapse, incl, epoch_type) :
def get_orbitxyz(orbelem) :
    ax, porb, ecc, epoch, Omega, periapse, incl, epoch_type = orbelem

    mjd0_JD = 2400000.5 
    if mjd0_JD < epoch :
        epoch_mjd = epoch - mjd0_JD
    else :
        epoch_mjd = epoch 

    ko = KeplerOrbit()
    ko.setOrbitalElements(ax, porb, ecc, epoch_mjd, Omega, periapse, incl);
    if epoch_type==0 :
        ko.setEpochTypeP()
    else :
        ko.setEpochTypeT()

    mjd_start = epoch_mjd
    mjd_end   = epoch_mjd + porb + 1e-6
    mjd_step  = porb/100.


    vmjd = np.arange(mjd_start, mjd_end, mjd_step)
    nv = len(vmjd)
    vx = np.zeros(nv)
    vy = np.zeros(nv)
    vz = np.zeros(nv)

    for idx in range(nv) :
        mjd = vmjd[idx]
        ko.calcPos(mjd)
        vx[idx] = ko.getPx()
        vy[idx] = ko.getPy()
        vz[idx] = ko.getPz()

    return vx, vy, vz
        

if __name__=="__main__" :

    #color_list = [
    #    ]
    #line
    linestyle_list = [
        "solid", "dotted", "dashed",
        ]

    parfname_list = [
        #"orbelems/herx1.param",
        #
        #"orbelems/cenx3.param",
        #"orbelems/velax1.param",
        #"orbelems/oao1657.param",
        #"orbelems/gx301m2.param",
        #
        "orbelems/rxj0520.param",
        "orbelems/4u0115.param",
        "orbelems/2s1553.param",
        "orbelems/v0332.param",
        "orbelems/ks1947.param",
        "orbelems/2s1417.param",
        "orbelems/exo2030.param",
        "orbelems/gs0834.param",
        "orbelems/a0535.param",
        "orbelems/gx304m1.param",
        "orbelems/xtej1946.param",
        "orbelems/groj1008.param",
        #
        ]

    Omega = 0
    incl  = 90.

    #mpl.style.use('default')
    #fig = plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    #plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    plt.subplots_adjust(left=0.1, right=0.72, top=0.9, bottom=0.1)

    idx = 0
    for parfname in parfname_list :
        if not os.path.isfile(parfname) :
            continue
        par = readfile_OrbitalElements(parfname)
        orbelem = (par[8], par[4], par[11], par[7], Omega, par[9], incl, par[6])
        target_name = par[0]
        print target_name, orbelem
        vx, vy, vz = get_orbitxyz(orbelem)
        
        #color_name = 'C%d'%(idx%10)
        linestyle_name = linestyle_list[idx/7]
        
        #plt.plot(vx, vy, ls=linestyle_name, label=target_name)
        plt.plot(vx, vz, ls=linestyle_name, label=target_name)
        idx+=1

    plt.xlim(-800, 800)
    plt.ylim(-800, 800)
    #plt.axes().set_aspect('equal', 'datalim')
    plt.axes().set_aspect('equal')
    plt.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0., fontsize=12)
    #plt.legend()
    
    #plt.savefig("orbit.png")
    #plt.savefig("orbit.eps")
    plt.savefig("orbit.pdf")
    plt.show()
