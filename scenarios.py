# -*- coding: utf-8 -*-
"""
Worst case functions

Copyright 2016 Cameron Dowd

This file is part of PropPy.

PropPy is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PropPy is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PropPy.  If not, see <http://www.gnu.org/licenses/>.

@author: c.dowd
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import data_dealings as data_dealings
import numpy.polynomial.polynomial as poly
from operator import itemgetter
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
import BLDC_model

import data_extraction

def static_max_current(C,atmosphere,battery,motor,propeller):
    """
    Finds maximum expected current for motor and propeller combination if held
    stationairy. Returns MOTOR current, not the current taken by the controller.
    """
    V_in=C*battery['V']
    
    rho=atmosphere['rho']
    D=propeller['diameter']  
    
    C_Q_interp, C_T_interp, bounds=data_dealings.prop_data_funcs(propeller)
    #Motor models
    Kv0=motor['Kv']*2*np.pi/60
    maxomega=battery['V']*Kv0 #Max omega at full voltage, won't ever be reached.
    
    Vinterp_func=None

    I_func=None
    
    #Initial bounds to search for the rotational speed, should contain decent guesses for all possible flight RPMs
    omega_search_min=0 #Don't bother with speeds that are highly unlikely to provide thrust, ie; outside bounds of surface
    
    omega_range=np.linspace(omega_search_min,maxomega,50)


    
    Q_range=data_dealings.prop_data_static_Q(omega_range,propeller,atmosphere)
    #Get motor voltage and current for the propeller curve
    if Vinterp_func is not None: 
        V_mot=Vinterp_func(Q_range,omega_range)
#        if V_mot!=V_mot: #If the values requested are outside the motor surface bounds, skip it
#            print(V_mot)
#            continue
        I_mot=I_func(Q_range) #Remember this is motor current, not system current, don't get too excited
    else:
        V_mot, I_mot, Temp= BLDC_model.Vsurf(Q_range,omega_range, motor, atmosphere)
    
    
    V_interp=interp1d(V_mot, omega_range, kind='quadratic')
    omega=[]
    #Find the omega that corresponds to each command voltage
    try:
        for V in V_in:
            try:
                omega.append(V_interp(V))
            except ValueError:
                omega.append(float('nan'))
    except TypeError:
        try:
            omega=V_interp(V_in)
        except ValueError:
            omega=float('nan')

    omega=np.array(omega)
    #Get torque and thrust for each of these real omega values
    Q=data_dealings.prop_data_static_Q(omega,propeller,atmosphere)
    thrust=data_dealings.prop_data_static_T(omega,propeller,atmosphere)

    
    if Vinterp_func is not None: 
        V_mot=Vinterp_func(Q,omega)
#        if V_mot!=V_mot: #If the values requested are outside the motor surface bounds, skip it
#            print(V_mot)
#            continue
        I_mot=I_func(Q) #Remember this is motor current, not system current, don't get too excited
    else:
        V_mot, I_mot, Temp= BLDC_model.Vsurf(Q,omega, motor, atmosphere)    
    
    rpm=omega*60/(2*np.pi)
    I_bat=I_mot*C #Get battery current

    return I_bat, rpm, thrust#, Q
    
if __name__ == "__main__":
    print("Run the main script")

    