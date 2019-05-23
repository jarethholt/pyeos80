#!/usr/bin/env python3
"""Equation of state of seawater EOS-80.

This library includes functions from EOS-80, a formulation of the equation of state of seawater. The primary purpose of the equation of state is to relate the enthalpy of seawater (and related thermodynamic quantities) to salinity, temperature, and pressure.
"""

# Import statements


# Constants
T68SCL = 1.00024  # Scale factor between ITS-90 and ITS-68 temperatures
DB2KPA = 1e-1  # Conversion factor from decibar to kilopascal

SECKPURE0 = (19652.21, 148.4206, -2.327105, 1.360477e-2, -5.155288e-5)
SECKPURE1 = (3.239908, 1.43713e-3, 1.16092e-4, -5.77905e-7)
SECKPURE2 = (8.50935e-5, -6.12293e-6, 5.2787e-8)
SECK02 = (54.6746, -0.603459, 1.09987e-2, -6.1670e-5)
SECK03 = (7.944e-2, 1.6483e-2, -5.3009e-4)
SECK12 = (2.2838e-3, -1.0981e-5, -1.6078e-6)
SECK13 = 1.91075e-4
SECK2 = (-9.9348e-7, 2.0816e-8, 9.1697e-10)

SMOW = (999.842594, 6.793952e-2, -9.095290e-3, 1.001685e-4, -1.120083e-6, 6.536332e-9)

DENS02 = (8.24493e-1, -4.0899e-3, 7.6438e-5, -8.2467e-7, 5.3875e-9)
DENS03 = (-5.72466e-3, 1.0227e-4, -1.6546e-6)
DENS04 = 4.8314e-4


# Helper functions
def poly1d(x,coefs):
    """Evaluate a polynomial in one variable.
    """
    y = 0.
    for coef in coefs[-1::-1]:
        y = y*x + coef
    return y


# Secant bulk modulus functions
def calseckpure0(t68):
    """Calculate the secant bulk modulus of pure water zeroth pressure term.
    """
    kw = poly1d(t68,SECKPURE0)
    return kw

def calseckpure1(t68):
    """Calculate the secant bulk modulus of pure water first pressure term.
    """
    aw = poly1d(t68,SECKPURE1)
    return aw

def calseckpure2(t68):
    """Calculate the secant bulk modulus of pure water second pressure term.
    """
    bw = poly1d(t68,SECKPURE2)
    return bw

def calseck0(spsu,t68):
    """Calculate the secant bulk modulus zeroth pressure term.
    """
    k0 = calseckpure0(t68)
    k2 = poly1d(t68,SECK02)
    k3 = poly1d(t68,SECK03)
    k = k0 + k2*spsu + k3*spsu**(1.5)
    return k

def calseck1(spsu,t68):
    """Calculate the secant bulk modulus first pressure term.
    """
    a0 = calseckpure1(t68)
    a2 = poly1d(t68,SECK12)
    a = a0 + a2*spsu + SECK13*spsu**(1.5)
    return a

def calseck2(spsu,t68):
    """Calculate the secant bulk modulus second pressure term.
    """
    b0 = calseckpure2(t68)
    b2 = poly1d(t68,SECK2)
    b = b0 + b2*spsu
    return b

def calseck(spsu,t90,pdb):
    """Calculate the secant bulk modulus of seawater.
    """
    # Convert temperatures and pressures
    t68 = t90 * T68SCL
    pkpa = pdb * DB2KPA
    
    # Get coefficients for the pressure polynomial
    k0 = calseck0(spsu,t68)
    a = calseck1(spsu,t68)
    b = calseck2(spsu,t68)
    coefs = (k0, a, b)
    
    # Calculate final bulk modulus
    k = poly1d(pkpa,coefs)
    return k


# Reference pure water (standard mean ocean water) functions
def caldsmow(t90):
    """Calculate density of pure water.
    """
    t68 = t90 * T68SCL
    dsmow = poly1d(t68,SMOW)
    return dsmow


# Atmospheric pressure functions
def caldens0(spsu,t90):
    """Calculate density of seawater at atmospheric pressure.
    """
    # Calculate normalized variables
    sqrts = spsu**.5
    t68 = t90 * T68SCL
    
    # Calculate coefficients of the salinity polynomial
    dsmow = caldsmow(t90)
    b = poly1d(t68,DENS02)
    c = poly1d(t68,DENS03)
    coefs = (dsmow, 0., b, c, DENS04)
    
    # Calculate final density
    dens0 = poly1d(sqrts,coefs)
    return dens0


# Density functions
def caldens(spsu,t90,pdb):
    """Calculate the density of seawater.
    """
    dens0 = caldens0(spsu,t90)
    seck = calseck(spsu,t90,pdb)
    pkpa = pdb * DB2KPA
    dens = dens0 / (1 - pkpa/seck)
    return dens
