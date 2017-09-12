#ifndef _keplerOrbit_H_
#define _keplerOrbit_H_

#include "atFunctions.h"


class KeplerOrbit {
 private :
  double m_a;     // semi-major axis (l/c sec)    
  double m_per;	  // orbital period (day)
  double m_ecc;	  // eccentricity		     
  double m_tau;	  // time of periastron passage (MJD)
  double m_Omega; // longitude of ascending node (rad)
  double m_w; 	  // argument of periastron (rad)
  double m_incl;  // orbit inclination (rad)
  int    m_epochtype; // P=Periaston epoch, T=epoch of l=M+omega=Pi/2  

  double m_nu;

  AtVect m_vel, m_pos;

  static const int EPOCH_TYPE_P = 0;
  static const int EPOCH_TYPE_T = 1;

 public:
  KeplerOrbit();
  virtual ~KeplerOrbit();
  
  int setOrbitalElements(double a, double per, double ecc, double tau, double Omega, double w, double incl);
  int setEpochTypeP();
  int setEpochTypeT();
  int printOrbitalElements();
  double getEccentricAnomaly(double mjd);
  double meanAnomaly(double mjd);
  double radius(double t, double E_in);

  double getDopplerFactor(double mjd);
  int calcPos(double mjd);
  int calcVel(double mjd);

  double getVx();
  double getVy();
  double getVz();

  double getPx();
  double getPy();
  double getPz();

  double getPorb();
  double getEpoch();
  double getPeriEpoch();

  double getAxsini();
  double getEccentricity();
  double getPeriArg();
  int getEpochType();
  
  //int xyzPos(double mjd, AtVect x);
  //int xyzVel(double mjd, AtVect vel);
};

#endif  /* _keplerOrbit_H_ */
