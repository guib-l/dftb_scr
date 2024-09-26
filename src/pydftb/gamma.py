#!/usr/bin/env python3
import __future__

# Import standard de pythons3
import os
import sys
import math
import copy

import numpy as np

###############################################################################
def Gamma(rab,Ua,Ub):
    """Gamma"""
    gammaE=lambda r,a,b:np.exp(-a*r)*( (0.5*b**4*a)*((a**2-b**2)**-2)-\
                    (b**6-3*b**4*a**2)*(r**1*(a**2-b**2)**3)**-1 )
    lim=lambda ua,ub:True if np.abs(ua-ub)<1e-6 else False
    tauA,tauB=3.2*Ua,3.2*Ub
    if rab<1e-6:
        if lim(Ua,Ub):
            return -0.5*(Ua+Ub)
        else:
            return -0.5*((tauA*tauB)*(tauA+tauB)**-1+\
                                    (tauA*tauB)**2*(tauA+tauB)**-3)
    elif lim(Ua,Ub):
        tauM=0.5*(tauA+tauB)
        return np.exp(-tauM*rab)*\
            (rab**-1+0.6875*tauM+0.1875*rab*tauM**2+0.0208333333333*rab**2*tauM**3)
    else:
        return gammaE(rab,tauA,tauB)+gammaE(rab,tauB,tauA)
    return ValueError



###############################################################################
def GammaPrime(rab,Ua,Ub):
    """GammaPrime"""
    gammaE=lambda r,a,b:-a*np.exp(-a*r)*( (0.5*b**4*a)*((a**2-b**2)**-2)-\
                    (b**6-3*b**4*a**2)*(r**1*(a**2-b**2)**3)**-1 ) + \
                    np.exp(-a*r)*( (b**6-3*a**4*b**2) * (r**2*(a**2-b**2)**3)**-1 )
    lim=lambda ua,ub:True if np.abs(ua-ub)<1e-6 else False
    
    if rab<1e-6:
        return 0.0
    elif lim(Ua,Ub):
        tauA,tauB=3.2*Ua,3.2*Ub
        tauM=0.5*(tauA+tauB)
        return -tauM*np.exp(-tauM*rab)*\
            (rab**-1+0.6875*tauM+0.1875*rab*tauM**2+0.0208333333333*rab**2*tauM**3) +\
            np.exp(-tauM*rab)*(-rab**-2+0.1875*tauM**2+2*0.0208333333333*rab*tauM**3) 
    else:
        tauA,tauB=3.2*Ua,3.2*Ub
        return gammaE(rab,tauA,tauB)+gammaE(rab,tauB,tauA)
    return ValueError


    