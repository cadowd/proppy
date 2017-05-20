# -*- coding: utf-8 -*-
"""
Functions to model BLDC electric motor from known motor constants.
The function allows a very complicated thermal model, allowing for the change
in constants with increasing motor temperature etc, however the values required
for the thermal model are hard to know. The model works close enough
if these extra motor values are set to 0.

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
import matplotlib.pyplot as plt
import csv
from itertools import islice
import data_dealings as data_dealings


def Vsurf(Qout,omega, motor, atmosphere): 
    """
    Function to predict the motor DC equivalent voltage and current as well
    as its running temperature for a given torque and motor speed. The inputs
    can be given as numpy arrays to speed the calculations.
    """
    Qout=np.array(abs(Qout))
    omega=np.array(omega)
    i0=motor['I0']
    T0=motor['Tref']
    R0=motor['R'] #Internal resistance at reference conditions in Ohms
    Kv0=motor['Kv']*2*np.pi/60 #Motor speed constant in rads/V, inverse of Kemf theoretically the same as KT for trapezoidal back EMF but is different due to non ideal back EMF waveform
    #KT0=1/Kv0*2/np.sqrt(3)
    KT0=1/Kv0 #Motor torque constant in Nm/A
    Ilimit=motor['Imax']
    alphar=motor['alphar'] #Resistivity temperature coefficient for copper windings
    alpham=motor['alpham'] #Magnet temperature coefficient for NdFeB magnet
    k1=motor['k1']
    k2=motor['k2']
    k3=motor['k3']
    Rthwm=motor['Rthwm'] #Thermal resistance winding to rotor (and magnets) [K/W],
    Rthma=motor['Rthma'] #Thermal resistance rotor to ambient [K/W],
    Rthws=motor['Rthws'] #Thermal resistance winding to stator [K/W], 
    Rthsa=motor['Rthsa'] #Thermal resistance stator to ambient [K/W], 
    
    if motor['use_thermal']:
            
        Rth1=Rthwm+Rthma
        Rth2=Rthws
        Rth3=Rthsa
    else:
        Rth1=0
        Rth2=0
        Rth3=0
#    print(Rthwm)

    #Initialise values
    Ta= atmosphere['Ta'] #ambient temperature
    Tw=Ta #winding temperature
    Ts=Ta #initial housing temperature
    Tm=Ta
    KT=KT0 #initial torque coefficient
    R=R0 #initial resistance
    n=0 #no iterations
    Pout=Qout*omega #Total mechanical power out
    
    Fa=omega/(10*Kv0) #assume no load current is measured at 10V, linear relationship
    
    Ifr=k3*Fa+k2*omega+k1*omega**2 #friction current (extra current required to overcome friction)
    I=Qout/KT+Ifr
    nmax=100 #max iterations
    maxTemp=2000 #max allowable motor temperature
    tol=10**-6 #Tolerance
    V=omega/Kv0+I*R #initial voltage
    """
    Here we iterate the temperature, updating the constants until the
    temperature stabilises. If the motor thermal constants are 0, only
    one iteration occurs.
    """
    #Iterate until convergence
    while True:
        n=n+1 #iteration number
        
        I=Qout/KT+Ifr #Current required to get operating torque including no load current

        Pin=V*I #Total power into motor

        Peldiss=I**2*R #Electrically dissipated power
                
        Pfr=Pin-Peldiss-Pout #The friction power is everything that isn't electrically dissipated
        
        if n>nmax: #Stop if it doesn't converge
            break
            print('Unconverged')
        
#        Tw=Peldiss*(Rthwh+Rthh0)+Pfr*Rthh0+Ta #Winding temperature, also stator temperature
        Ts1=Ts
        if not (Rth1==0 and Rth2==0 and Rth3==0):
            Tw=(Ta*(Rth1+Rth2+Rth3)+Rth1*(Pfr*Rth3+Peldiss*(Rth2+Rth3)))/(Rth1+Rth2+Rth3)  #Winding temperature, also stator temperature
            Ts= (Ta*Rth2+Tw*Rth3+Rth2*Rth3*Pfr)/(Rth3+Rth2)#Stator temperature
            Tm=Rthma*(Tw+Ta)/Rth1+Ta #magnet temperature, also rotor temperature
        try:
            if max(Tw)>maxTemp:
                print('Over temperature limit')
                break
        except TypeError:
            if Tw>maxTemp:
                break
        
        deltaTw=Tw-T0 #Increase in winding temperature (FROM NOMINAL VALUE, NOT AMBIENT)
        deltaTmag=Tm-T0 #increase in magnet temperature
    
        
        R=R0*(1+deltaTw*alphar) #Resistance at temperature
        Kv=Kv0*(1+deltaTmag*alpham) #Speed constant at new temperature
        KT=1/Kv
        
        V=omega/Kv+I*R #Voltage
        
        
        try:
            if max(abs(Ts1-Ts))<tol:
                break
        except TypeError:
            if abs(Ts1-Ts)<tol:
                break 
        
#    print(Ts1)
#    print(Tm)
    return(V,I,Ts1)

def opt_kv(omega,Q, Kv):
    I0=3.0
    R=0.115
    Kv=Kv/60*2*np.pi
    Pel=Q*omega+omega*I0/Kv + R*Q**2*Kv**2 + 2*R*Kv*I0 + R*I0**2
    
    return Pel

if __name__ == "__main__":
    print("Run the other one")
    
