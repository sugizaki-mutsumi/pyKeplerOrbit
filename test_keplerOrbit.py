#!/usr/bin/env python

from keplerOrbit import KeplerOrbit


#name         = A 0535+26   // Object Name
ra           = 84.727383  #// (deg)
dec          = 26.315784  #// (deg)
#binary       = Y   // Y=yes, N=no
pbinary      = 111.100000  #// Orbital period (days)
pbdot        = 0.000000    #// Derivative of Orbital period (days/days)
#epoch_type   = P   // P=Periaston epoch, T=epoch of l=M+omega=Pi/2
binaryepoch  = 2453613.500000  #// Epoch (JD)
axsini       = 267.000000  #// Projected semimajor axis (light-sec)
periapse     = 130.000000  #// Argument of periapse (deg)
apsidalrate  = 0.000000  #// Rate of change of arg. of periapse (deg/year)
eccentricity = 0.470000  #// Eccentricity
#egress       = 0.000000  // Binary phase at end of eclipse (cycles)
#ingress      = 1.000000  // Binary phase at begining of eclipse (cycle)

epoch_mjd = binaryepoch - 2400000.5

Omega=0.0
incl =90.0
ko = KeplerOrbit()
#ko.setOrbitalElements(axsini, pbinary, eccentricity, epoch_mjd, Omega, periapse,incl);
ko.init_OrbitalElements(axsini, pbinary, eccentricity, epoch_mjd, Omega, periapse,incl);

mjd_start = 55050
mjd_end   = 56000
mjd_step  = 1

mjd = mjd_start
while mjd<mjd_end :
    dopplerfact = ko.getDopplerFactor(mjd)
    print(mjd, dopplerfact)
    mjd+=mjd_step

