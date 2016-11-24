# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 12:04:48 2016
Flying wing drag polar estimator for low Re model scales
@author: c.dowd
"""

import numpy as np
import consumption_functions
import matplotlib.pyplot as plt

def getSpans(U_stall, plane, atmosphere):
    """
    Returns the drag of a given plane in a given atmosphere at the specified
    flight speed. The plane is modelled through a quadratic approximation
    to its drag coefficient to lift coefficient polar.
    """
    mass=payload/payload_frac
    aero_factor=plane['aero_factor']
    rho=atmosphere['rho']
#    mass=plane['mass']
    
    
    
    L=9.8*mass
    tc_ratio=0.2
    
    
    
    CL_max=1.0+aero_factor #range from 1 to 2 (best possible with flaps)
    
    CL_max_swept=CL_max*np.cos(theta_sweep)
    
    U_thrown_min=5
    
    static_acc_min=5
    
    allowable_drop=0.75
    
    drop_time=np.sqrt(2*allowable_drop/9.81)
    U_stall=U_thrown_min+static_acc_min*drop_time
    print(plane_mass)
    
    S_ref=2*L/(rho*U_stall**2*CL_max_swept)
    print(S_ref)
    
    c_root=body_length
    
    AR_max=4*S_ref/c_root**2 #50.0

    AR_min=S_ref/c_root**2 #2
    b_min=np.sqrt(S_ref*AR_min)
    b_max=np.sqrt(S_ref*AR_max)
    
#    print(b_max/AR_max)
    AR=AR_min + AR_factor*(AR_max-AR_min)
    b=np.sqrt(S_ref*AR)
    c_root=body_length
    c_tip=2*S_ref/b-c_root

def dragFunc(U, plane, atmosphere):
    """
    Returns the drag of a given plane in a given atmosphere at the specified
    flight speed. The plane is modelled through a quadratic approximation
    to its drag coefficient to lift coefficient polar.
    """
    
    
    print(b)
    

        
    Z=2*np.cos(theta_sweep)
    K=1+Z*(tc_ratio)+100*(tc_ratio)**4 #Form factor of wing, Shevell 1989
    
    print(K)
    Re=U*np.sqrt(S_ref)*rho/atmosphere['mu'] #Reynolds number
    C_f=0.455/(np.log10(Re)**2.58) #Schlichting Prandtl skin friction    
#    print(C_f)
    S_wet=S_ref*2*1.02    
    
    Cd_p=Q/Sref*K*C_f*S_wet
    print(Cd_p)
    #Add for fuselage and vertical surfaces
    
    lambda_t=c_tip/c_root - 0.357 + 0.45*np.exp(0.0375*theta_sweep) #Taper ratio with sweep correction
    func_lambda=0.0524*lambda_t**4 - 0.15*lambda_t**3 + 0.1659*lambda_t**2 - 0.0706*lambda_t + 0.0119
    e_th=1/(1+func_lambda*AR)
    
    #Correction factors for the Oswald form factor
    ke_D0=0.9 #Rough estimation for general aviation
    ke_WL=1.00 #Should be accounted for properly maybe
    
    ke_f= 1-2*(body_height/b) #Accounts for fuselage
    
    e_total=e_th*ke_f*ke_D0*ke_WL #total Oswald form factor
    print(e_total)
    
    Cl=2*L/(rho*U**2*S_ref)
    
    Cd_if=1/(np.pi*e_total*AR)
    print(Cd_if)
    
    C_d=Cd_p + Cd_if*Cl**2
#    C_d=0.0174 + 0.1156*Cl**2

    
#    Re=U*np.sqrt(plane['ref_area'])*atmosphere['rho']/atmosphere['mu'] #Reynolds number
    drag=C_d*0.5*rho*plane['ref_area']*U**2
    return drag
    

if __name__ == "__main__":
    L=0.7 #Wing length
    c=0.2 #Wing chord length
    Sref=0.3259#+0.2
    
    b_t= 1.2 #total span of craft
    
    Q=1.2 #to 1.3 depends on other protruding parts
    
    payload=0.15 #payload mass
    theta_sweep=30*np.pi/180
    body_length=0.4
    body_height=0.07
#    rho=1.15
    
    payload_frac=0.125 #RESEARCH
    
    
    aero_factor= 1.0 #0 to 1
    AR_factor=0.2 #0 to 1
    
    plane1={
    'C_d':0.04, 
    'ref_area':0.3259,
    'mass': 1.15, #mass in kilos
    'C_1': 0.1156,
    'C_2': -0.0069,
    'C_3': 0.0174,
    'loss_factor': 1.0, #Assumed to scale with reynolds
    'loss_ref_length': 1,
    'theta_sweep': 30,
    'body_length': 0.4,
    'aero_factor':0.5
    }
    
    atmosphere={
    'rho':1.10, #
    'mu':1.8e-5, #dynamic viscosity
    'Ta': 25, #Ambient temperature
    }    
    
    plane_mass=payload/payload_frac    
    
    U_range=np.linspace(5,20,100)
    
    drag=dragFunc(U_range,plane1,atmosphere)
    
    drag_real=consumption_functions.dragFunc(U_range,plane1,atmosphere)
    
    plt.plot(U_range, drag, label='Estimate')
    
    plt.plot(U_range, drag_real, label='real')
    plt.legend(loc=4)
    plt.show()
