#!/usr/bin/env python
# Christoph Weniger, 6 Dec 2012
#
# Numbers loosely oriented on HESS-II (needs update)
# Information on ROI must be updated as well (J-value, size)

from __future__ import division
from numpy import *
from scipy import integrate, stats

# Proton background (for references see 1106.1874)
def bg_protons(x):
# astro-ph/0210453
  return 8.73e-9*(x/1000)**-2.71 # 1/cm^2 s GeV sr
e_rec = 1e-2 # proton rejection factor adopted for HESS (same as in 1207.6773)
shift = 3 # proton energy shift (due to reconstruction as gamma rays)

# Exposure etc at 130 GeV
exposure =  .01 * 1e5**2 # .01 km^2 in cm^2
time = 5*3600 # 50 hours

# ROI properties
box_size = 4.0 # length and height in deg
sr_roi = radians(box_size)**2 # ROI size in sr
#Jvalue = 2.24664205472e+22 # l.o.s. integral over rho^2, integral over d\Omega [GeV^2 cm^-5 sr]
Jvalue = 5.20952217306e+22 # l.o.s. integral over rho^2, integral over d\Omega [GeV^2 cm^-5 sr]

# Relevant energies and energy resolution
Erange = 110, 150 # Energy range
Eline = 130 # Line energy
Eres = 0.2 # Energy resolution (68% containment)

# Electron background (cannot be rejected efficiently)
# Abdo et al., PRL 102 (2009) 181101
def bg_electrons(x):
  return 1.5e-5*(x/10.)**-3.0 # 1/cm^2 s GeV sr


def Calc_J_Factor():
  import MC
  #profile = ('NFW',23.5,1.0)
  profile = ('EIN',20.0,0.17)
  MC.Gen_Annihilation_Map(box_size, 200, profile,'ein_test.pickle')

# Run to calc new J-Factor.  Takes a few minutes.
#Calc_J_Factor()


# Effective proton background (after energy shift and rejection)
def bg_protons_effective(x):
  return bg_protons(x*shift)*e_rec*shift

# This gives dN/dE integrated over the energy window, assuming annihilation
# into gamma gamma
def line_photons_within_ebin():
  return 2*integrate.quad(lambda x: stats.norm.pdf(x, Eline, Eline*Eres), Erange[0], Erange[1])[0]

def main():
  print 'Number of background photons:',
  dNbg = lambda x: (bg_electrons(x) + bg_protons_effective(x))*exposure*time*sr_roi
  dNbg_electrons = lambda x: (bg_electrons(x))*exposure*time*sr_roi
  Nbg = integrate.quad(dNbg, Erange[0], Erange[1])[0]
  Nbg_electrons = integrate.quad(dNbg_electrons, Erange[0], Erange[1])[0]
  print Nbg
  print 'Number of background electrons only:', Nbg_electrons
  
  

  print 'Number of signal photons:',
  Nsig = Jvalue*time*exposure/8/pi*1.27e-27/Eline**2*line_photons_within_ebin()
  print Nsig

  print "Signal to background:", Nsig/Nbg

  # For cross-checks:
  # print Nbg/exposure/time/(Erange[1]-Erange[0])*Eline**2
  # print bg_electrons(100)*100**3*1e4

if __name__ == '__main__':
  main()
