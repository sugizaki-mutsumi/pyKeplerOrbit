#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>

#include "atFunctions.h"

#include "keplerOrbit.h"

KeplerOrbit::KeplerOrbit(){
}

KeplerOrbit::~KeplerOrbit(){
}
  
int KeplerOrbit::setOrbitalElements(double a, double per, double ecc, double tau, double Omega, double w, double incl){
  m_a   = a;                    
  m_per = per;               
  m_ecc = ecc;               
  m_tau = tau;               
  m_Omega = Omega*DEG2RAD;   
  m_w     = w*DEG2RAD;         
  m_incl  = incl*DEG2RAD;     
  m_epochtype = EPOCH_TYPE_P; // default
  
  m_nu = 2.0*M_PI/m_per;     

  return 1;
}
    
int KeplerOrbit::printOrbitalElements(){
  std::cout << "axsini = " << m_a << std::endl;                    
  std::cout << "Orbital period = " << m_per    << std::endl;                    
  std::cout << "Eccentricity   = " << m_ecc    << std::endl;                    
  std::cout << "Epoch          = " << m_tau    << std::endl;                    
  std::cout << "Omega          = " << m_Omega*RAD2DEG   << std::endl;                    
  std::cout << "w              = " << m_w*RAD2DEG       << std::endl;                    
  std::cout << "Inclination    = " << m_incl*RAD2DEG    << std::endl;                    
  std::cout << "Epoch type     = " << m_epochtype  << std::endl;                    

  return 1;
}
    
int KeplerOrbit::setEpochTypeT(){
  m_epochtype = EPOCH_TYPE_T; /// default
  return 1;
}
int KeplerOrbit::setEpochTypeP(){
  m_epochtype = EPOCH_TYPE_P; /// default
  return 1;
}


double KeplerOrbit::getEccentricAnomaly(double mjd){
  double M, E;
  M = meanAnomaly(mjd);
  //E = getE(M, m_ecc);
  atKepler(M, m_ecc, &E);
  return E;
}

double KeplerOrbit::meanAnomaly(double mjd){
  //return m_nu*(t-m_tau); /// bug in ver 2015.09
  if(m_epochtype==EPOCH_TYPE_T){
    return m_nu*(mjd-m_tau)-m_w+M_PI/2.0; 
  } else {
    return m_nu*(mjd-m_tau); /// normal periastron epoch
  }
};
  
double KeplerOrbit::radius(double mjd, double E_in){
  double E;
  if(E_in<10.0) E = getEccentricAnomaly(mjd);
  else E = E_in;
  return m_a * (1.-m_ecc*cos(E));
}

//int KeplerOrbit::xyzPos(double mjd, AtVect x){
int KeplerOrbit::calcPos(double mjd){
  double E, r, f, wf;
  E = getEccentricAnomaly(mjd);
  r = radius(mjd, E);
  ////f = arctan( sqrt(1.-m_ecc*m_ecc) * tan(E/2.) ) * 2.0;
  //f = atan( sqrt(1.-m_ecc*m_ecc) * tan(E/2.) ) * 2.0; // Bug
  f = atan( sqrt((1.+m_ecc)/(1.-m_ecc)) * tan(E/2.) ) * 2.0; // Bug fix
  wf = m_w + f;
    
  m_pos[0] = r*( cos(m_Omega)*cos(wf) - sin(m_Omega)*sin(wf)*cos(m_incl) );
  m_pos[1] = r*( sin(m_Omega)*cos(wf) + cos(m_Omega)*sin(wf)*cos(m_incl) );
  m_pos[2] = r*( sin(wf)*sin(m_incl) );
  
  return 1;
}

//int KeplerOrbit::xyzVel(double mjd, AtVect vel){
int KeplerOrbit::calcVel(double mjd){
  double E;
  double l1, l2, m1, m2, n1, n2, b, r, nar;

  E = getEccentricAnomaly(mjd);
  l1 = cos(m_Omega)*cos(m_w) - sin(m_Omega)*sin(m_w)*cos(m_incl);
  l2 = -cos(m_Omega)*sin(m_w) - sin(m_Omega)*cos(m_w)*cos(m_incl);
  m1 = sin(m_Omega)*cos(m_w) + cos(m_Omega)*sin(m_w)*cos(m_incl);
  m2 = -sin(m_Omega)*sin(m_w) + cos(m_Omega)*cos(m_w)*cos(m_incl);
  n1 = sin(m_w)*sin(m_incl);
  n2 = cos(m_w)*sin(m_incl);
  b = m_a * sqrt(1. - m_ecc*m_ecc);
  r = radius(mjd, E);
  nar = (m_nu/86400.0) * m_a / r;

  m_vel[0] = nar * ( b*l2*cos(E) - m_a*l1*sin(E) );
  m_vel[1] = nar * ( b*m2*cos(E) - m_a*m1*sin(E) );
  m_vel[2] = nar * ( b*n2*cos(E) - m_a*n1*sin(E) );

  return 1;
};

double KeplerOrbit::getDopplerFactor(double mjd){
  double dopplerfact;
  calcVel(mjd);
  dopplerfact = sqrt( (1+m_vel[2])/(1-m_vel[2]) );
  return dopplerfact;
}

double KeplerOrbit::getVx(){return m_vel[0];};
double KeplerOrbit::getVy(){return m_vel[1];};
double KeplerOrbit::getVz(){return m_vel[2];};

double KeplerOrbit::getPx(){return m_pos[0];};
double KeplerOrbit::getPy(){return m_pos[1];};
double KeplerOrbit::getPz(){return m_pos[2];};

double KeplerOrbit::getPorb(){return m_per;};
double KeplerOrbit::getEpoch(){return m_tau;};
double KeplerOrbit::getPeriEpoch(){
  if(m_epochtype==EPOCH_TYPE_T){
    return m_tau + m_per*(m_w/(2*M_PI)-0.25);
  } else {
    return m_tau; /// normal periastron epoch
  }
};

double KeplerOrbit::getAxsini(){return m_a;};
double KeplerOrbit::getEccentricity(){return m_ecc;};
double KeplerOrbit::getPeriArg(){return m_w;};
int    KeplerOrbit::getEpochType(){return m_epochtype;};
