# -*- coding: utf-8 -*-
"""
Created on Sat May 20 20:37:10 2017

@author: lis-15-15
"""
import numpy as np

def CT_model(prop_ratio):
    """
    Estimates thrust coefficient of a propeller based on its diameter/pitch
    ratio. The model is based off a linear fit to the entirety of APC propeller
    data from the University of Illinois propeller database.
    """
    model_coeffs=[-0.01733052, 0.14146734]
    CT_fit_func=np.poly1d(model_coeffs)
    CT=CT_fit_func(prop_ratio)
    if CT<0:
        CT=0.0001
    return CT

def CQ_model(prop_ratio):
    """
    Estimates torque coefficient of a propeller based on its diameter/pitch
    ratio. The model is based off an exponential fit to the entirety of APC propeller
    data from the University of Illinois propeller database.
    """
    model_coeffs=[ 0.06778287,  1.83984961,  0.00490218]
    CQ_fitted=model_coeffs[0]*np.exp(-model_coeffs[1]*prop_ratio)+model_coeffs[2]
    
    return CQ_fitted

def slip_factor_model(prop_ratio):
    """
    Estimates torque coefficient of a propeller based on its diameter/pitch
    ratio. The model is based off an exponential fit to the entirety of APC propeller
    data from the University of Illinois propeller database.
    """
    model_coeffs=[ 1.76372577,  4.15848726,  1.15502056]
    slip_factor_fitted=model_coeffs[0]*np.exp(-model_coeffs[1]*prop_ratio)+model_coeffs[2]
    
    return slip_factor_fitted

def opt_kv(V,omega_max,Q_max):
    R=0.2 #Rough guess, should be refined based on data
    if 4*Q_max*R*omega_max>V**2:
        return "Error"
    Kv_opt=(V-np.sqrt(V**2-4*Q_max*R*omega_max))/(2*Q_max*R)*30/np.pi
    
    return Kv_opt