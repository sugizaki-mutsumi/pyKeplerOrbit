#ifndef _keplerOrbit_H_
#define _keplerOrbit_H_

//#include "atFunctions.h"
static const double DEG2RAD=0.017453292519943295769;
static const double RAD2DEG=57.295779513082320877;
typedef double AtVect[3];

// for Kepler solver
static const double EPS = 1e-15;
static const int IMAX = 50;


class KeplerOrbit {
 private :
  double m_a;     // semi-major axis (l/c sec)    
  double m_per;	  // orbital period (day)
  double m_pdot;  // derivative of orbital period (day/day)

  double m_ecc;	  // eccentricity		     
  double m_tau;	  // epoch in MJD : time of periastron (P) or l=M+omega=Pi/2 (T)
  double m_Omega; // longitude of ascending node (rad)
  double m_w; 	  // argument of periastron (rad)

  double m_arate; // Rate of change of arg. of periapse (deg/year)
  double m_incl;  // orbit inclination (rad)
  int    m_epochtype; // P=Periaston epoch, T=epoch of l=M+omega=Pi/2  

  double m_nu;

  AtVect m_vel, m_pos;

  static const int EPOCH_TYPE_P = 0;
  static const int EPOCH_TYPE_T = 1;

  
 public:
  KeplerOrbit();
  virtual ~KeplerOrbit();

  int init_OrbitalElements(double a, double per, double ecc, double tau, double Omega, double w, double incl);
  int init_OrbitalElements2(double a, double per, double pdot, double ecc, double tau, double Omega, double w, double arate, double incl);

  int setOrbitalPeriod(double per);
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
