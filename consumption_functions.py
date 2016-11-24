# -*- coding: utf-8 -*-
"""
Created on Nov 20 10:04:13 2016

Contains the functions for calculating the power consumption of a
given motor and propeller combination for a plane. For an explanation of roughly
how they work see the PropPy use notes. For a detailed explanation
of how they work read this code and clear your afternoon.

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
#import scenarios
import numpy.polynomial.polynomial as poly
from operator import itemgetter
import itertools
from matplotlib import _cntr as cntr
import scenarios
import BLDC_model

def funky_generator(atmosphere,battery,motor,propeller):
    
    #Load constant models. This loads all the data we have for this combination, and is the only time we access the files.
    #propeller models
    C_Q_interp, C_T_interp, bounds=data_dealings.prop_data_funcs(propeller)
    #Motor models
    Kv0=motor['Kv']*2*np.pi/60
    maxQ=100/Kv0 #Max torque at 100 amps. Thats an insane value but whatever.
    maxomega=25*Kv0 #Max omega at 25 volts. Quite a lot.
    Vinterp_func=None
    mot_bounds=[0,maxQ, 0, maxomega]
    I_func=None
    
    prop_funcs={'C_Q_interp':C_Q_interp,
                         'C_T_interp':C_T_interp,
                         'bounds':bounds,
                         'propeller':propeller
                         }
    motor_funcs={'Vinterp_func':Vinterp_func,
                         'I_func':I_func,
                         'bounds':mot_bounds,
                         'motor': motor
                         }
    return motor_funcs, prop_funcs

def P_el_calculator(motor_funcs,prop_funcs,plane, battery, atmosphere):
    """
    Returns electrical power consumption curve of given plane and motor and propeller functions (generated by
    the generation function). First it creates a flight envelope of all possible propeller rotation speeds
    and flight speeds that result in steady flight. This results in a list of torque values for steady flight.
    At each of these torque and rotational speeds the motor's consumption is calculated from the motor voltage
    surface. The smoothness and general likeability and character of the voltage surface depends strongly on the
    motor data collected, if you don't like the results this function probably doesn't like your data. Check
    that no extreme values show up in the csv files and that the columns values actually contain what they claim
    and aren't mixed up by mistake (as the test bench often changes port numbers on its own just to see if you're
    paying attention).
    """
    propeller=prop_funcs['propeller']
    C_T_interp=prop_funcs['C_T_interp']
    Vinterp_func=motor_funcs['Vinterp_func']
    I_func=motor_funcs['I_func']
    V_bat=battery['V']
    
    D=propeller['diameter']
    propeller_funcs=prop_funcs
    
    #Find the location of minimum drag to allow us to better resolve the drag bucket
    drag_min_vel=sp.optimize.minimize_scalar(lambda x : dragFunc(x, plane, atmosphere), bounds=(0.5, 100), method='bounded')
    drag_min_vel=drag_min_vel.x
    drag_min=dragFunc(drag_min_vel, plane, atmosphere)
#    print(drag_min_vel)
    
    #Make the bounds of our buckety buddy
    bucket_r=drag_min_vel+drag_min_vel*0.05
    bucket_l=drag_min_vel-drag_min_vel*0.2
    

    #Initial bounds to search for the rotational speed, should contain decent guesses for all possible flight RPMs
    omega_search_min=2/(prop_funcs['bounds'][3]*D)*2*np.pi #Don't bother with speeds that are highly unlikely to provide thrust, ie; outside bounds of surface
#    print(omega_search_min)
    omega_search_max=motor_funcs['bounds'][3]* V_bat/12#We also scale it by the battery voltage, as the surface max is based on a maximum of 25 volts
#    print(omega_search_max)
    if omega_search_max!=omega_search_max: #Ensure it exists
        print('Maximum command not enough to turn motor.')
        flight_status=False
        results_list=[{'omega':0,
                             'U':0,
                             'I_mot':0,
                             'P_el': 0,
                             'P_mech_prop':0,
                             'P_mech_prop':0,
                             'flight_status':flight_status}]
        return results_list
    
    U_search_min=7#minimum velocity we even consider
    U_search_max=prop_funcs['bounds'][3]*omega_search_max/np.pi*0.5*D  #Again, at some speed the plane will not produce thrust, we dont bother searching above this speed.
    #The max velocity is a little extreme but it is theoretically possible if the motor is fuckin ripped brah

    #Resolution of search curves    
    #To find intersections faster, searches are done using interpolated values across a range.
    search_res=50
    omega_range=np.linspace(omega_search_min,omega_search_max,search_res)
    bucket_range=np.linspace(bucket_l,bucket_r, search_res)
    U_range = np.linspace(U_search_min, U_search_max, search_res)
    U_range=np.sort(np.append(U_range, bucket_range))

    #Calculate the envelope of possible flight speeds and propeller rotations
    try:
        omega_flight, U_flight, Q_flight=flightEn(U_range, omega_range, propeller_funcs, propeller, plane, atmosphere)
    except TypeError:
        print("Prop/motor combination unable to sustain level flight.")
        results_list=[]
        return results_list
    #This is our 'big boy' interpolation, we resolve each point more finely to improve the result.
    #Now we find the motor consumption at this torque and rpm
    results_list=[]
    for i, omega in enumerate(omega_flight):
        Q=Q_flight[i]
        U=U_flight[i]
        if Vinterp_func is not None: 
            V_mot=Vinterp_func(Q,omega)
            if V_mot!=V_mot: #If the values requested are outside the motor surface bounds, skip it
                print(V_mot)
                continue
            I_mot=I_func(Q) #Remember this is motor current, not system current, don't get too excited
        else:
            V_mot, I_mot, Temp= BLDC_model.Vsurf(Q,omega, motor_funcs['motor'], atmosphere)
        P_el=V_mot*I_mot
        Th=propT(omega, U, C_T_interp, prop_funcs['bounds'], propeller,atmosphere)[0]
        P_mech_prop=Th*U
        P_mech_shaft=Q*omega
        if V_mot<V_bat: #Don't go above battery voltage, also ensures that the voltage is a real number and not something odd
            results_list.append({'omega':omega,
                                 'U':U,
                                 'I_mot':I_mot,
                                 'P_el': P_el,
                                 'P_mech_prop':P_mech_prop,
                                 'P_mech_shaft':P_mech_shaft})
        
#    print(results_list)
    return results_list#Pel, U_val, rpm, I_in, P_mech_shaft, P_mech_prop, status

def flightEn(U_range, omega_range, propeller_funcs, propeller, plane, atmosphere):
    """
    Finds all possible propeller rotation speeds, flight speeds and propeller
    torque values for a given propeller and plane within a range of speeds
    and rotation speeds.
    """
    C_T_interp=propeller_funcs['C_T_interp']
    C_Q_interp=propeller_funcs['C_Q_interp']
    bounds=propeller_funcs['bounds']
    #We make a grid of rotational speeds and flight speeds representing the possible flight envelope
    X, Y = np.meshgrid(omega_range, U_range)
    #We calculate the thrust surface at every single one of these shits
    Ts = np.array([propT(x,y, C_T_interp, bounds, propeller,atmosphere) for x,y in zip(np.ravel(X), np.ravel(Y))])
    Tsurf = Ts.reshape(X.shape)
    #Z=T#C_T_interp(X,Y)
    
    # Calculate the drag surface
    Ds= np.array([dragFunc(y, plane, atmosphere) for x,y in zip(np.ravel(X), np.ravel(Y))])
    Dsurf = Ds.reshape(X.shape)

    # Take the difference between the two surface heights and find the contour
    # where that surface is zero.
    diffSurf = Dsurf - Tsurf;

    #The contour represents the possible flight speeds and rotation speeds that result in level flight
    #This includes reynolds number and advance ratio effects on the thrust coefficient.
    c = cntr.Cntr(X, Y, diffSurf)
    # trace a contour at z == 0.0
    res = c.trace(0.0)
    
    # result is a list of arrays of vertices and path codes
    nseg = len(res) // 2
    if nseg==0:
        print('No intersection of drag and thrust surface for given U and omega range found.')
        return
    segments, codes = res[:nseg], res[nseg:]
    
    #Suck out the points on this contour
    omega_flight = segments[0][:,0]
    U_flight = segments[0][:,1]
    if np.isnan(omega_flight).any():
        print('UGH FUCK')
        print(omega_range)
        print(U_range)
        print(Dsurf)
    if np.isnan(U_flight).any():
        print('SHIT')
#    print('length is ' + str(len(omega_flight)))

    #Improve our results by repeating within the limits of the outputted flight envolope
    #We expand the envelope slightly to ensure we capture the edges properly.
    omega_range=np.linspace(0.9*min(omega_flight), max(omega_flight)*1.1, len(omega_range))
    U_range=np.linspace(0.9*min(U_flight), max(U_flight)*1.1, len(U_range))
    X, Y = np.meshgrid(omega_range, U_range)
    Ts = np.array([propT(x,y, C_T_interp, bounds, propeller,atmosphere) for x,y in zip(np.ravel(X), np.ravel(Y))])
    Tsurf = Ts.reshape(X.shape)
    Ds= np.array([dragFunc(y, plane, atmosphere) for x,y in zip(np.ravel(X), np.ravel(Y))])
    Dsurf = Ds.reshape(X.shape)

    diffSurf = Dsurf - Tsurf;
    c = cntr.Cntr(X, Y, diffSurf)
    res = c.trace(0.0)
    
    nseg = len(res) // 2
    if nseg==0:
        print('No intersection of drag and thrust surface for given U and omega range found.')
        return
    segments, codes = res[:nseg], res[nseg:]
    #Suck out the points on this contour
    omega_flight = segments[0][:,0]
    U_flight = segments[0][:,1]
    #For each point we find the torque
#    print('resolved length is ' + str(len(omega_flight)))
    Q_flight= np.array([propQ(x,y, C_Q_interp, bounds, propeller,atmosphere) for x,y in zip(omega_flight, U_flight)])

    
    return omega_flight, U_flight, Q_flight.reshape(U_flight.shape)

def dragFunc(U, plane, atmosphere):
    """
    Returns the drag of a given plane in a given atmosphere at the specified
    flight speed. The plane is modelled through a quadratic approximation
    to its drag coefficient to lift coefficient polar.
    """
    L=9.8*plane['mass']
    Cl=2*L/(atmosphere['rho']*U**2*plane['ref_area'])
    C_d=plane['C_1']*Cl**2+plane['C_2']*Cl+plane['C_3']
    drag=C_d*0.5*atmosphere['rho']*plane['ref_area']*U**2
    return drag

def find_zeros(diff_array):
    i_diffs=[]
    a_diffs=[]
    for i in range(len(diff_array) - 1):
        if diff_array[i] == 0. or diff_array[i] * diff_array[i + 1] < 0.:
            #Here we're looping over every value in a difference function, checking if it changes signs.
            # crossover at i
            i_diffs.append(i)
            a=abs(diff_array[i])/(abs(diff_array[i])+abs(diff_array[i+1]))
            a_diffs.append(a)
    return i_diffs, a_diffs

def propQ(omega,U, C_Q_interp_surface, bounds, propeller,atmosphere):
    """
    Returns list of torques for a propeller given a CQ interpolation surface and its bounds
    for a given omega and flight speed U.
    If the requested values are outside the bounds of the interpolation surface,
    the function returns the value of the surface at its limits.
    """
    Jrange=U/(omega/(2*np.pi)*propeller['diameter'])
    #Calculate blade reynolds number at these conditions
    tip_vel_rot=omega*propeller['diameter']/2
    tip_vel_tot=np.sqrt(tip_vel_rot**2+U**2)
    Rerange=tip_vel_tot*propeller['diameter']/(atmosphere['mu'])
    Q=[]
    if isinstance( Jrange, (float, int) ): #make J a list if we've fed in floats.
            Jrange=[Jrange]
            Rerange=[Rerange]
    for i,J in enumerate(Jrange):
        Re=Rerange[i]
        #Limit ourselves to values in bounds of surface, if its not on the surface i don't want to hear about it
        if J>bounds[3]:
            J=bounds[3]
        if Re>bounds[1]:
            Re=bounds[1]
        if J<bounds[2]:
            J=bounds[2]
        if Re<bounds[0]:
            Re=bounds[0]

        #Interpolate the torque coefficient using a 2D interpolation of the advance ratio and reynolds number
        C_Q=C_Q_interp_surface([J,Re])[0]
        if np.isnan(C_Q): #The interpolation returns nan if the initial points are outside range.
            print('ERROR REQUESTED VALUE OUTSIDE BOUNDS FOR Q')
#            print(omega)
#            print(U)
        #Check if omega is a value or array, so we can use this function in both cases.
        if not isinstance( omega, (float, int) ):
#            print(i)
#            print(omega)
            omega_val=omega[i]
        else:
            omega_val=omega
        Q.append(C_Q*atmosphere['rho']*(omega_val/(2*np.pi))**2*propeller['diameter']**5)
    return np.array(Q)
    
def propT(omega,U, C_T_interp_surface, bounds, propeller,atmosphere):
    """
    Returns list of thrusts for a propeller given a CT interpolation surface and its bounds
    for a given omega and flight speed U.
    If the requested values are outside the bounds of the interpolation surface,
    the function returns the value of the surface at its limits.
    """
    Jrange=U/(omega/(2*np.pi)*propeller['diameter'])
    #Calculate blade reynolds number at these conditions
    tip_vel_rot=omega*propeller['diameter']/2
    tip_vel_tot=np.sqrt(tip_vel_rot**2+U**2)
    Rerange=tip_vel_tot*propeller['diameter']/(atmosphere['mu'])
    T=[]
    if isinstance( Jrange, (float, int) ): #make J a list if we've fed in floats.
            Jrange=[Jrange]
            Rerange=[Rerange]
    for i,J in enumerate(Jrange):
#        print(J)
        Re=Rerange[i]
        if J>bounds[3]:
            J=bounds[3]
        if Re>bounds[1]:
            Re=bounds[1]
        #Interpolate the torque coefficient using a 2D interpolation of the advance ratio and reynolds number
        C_T=C_T_interp_surface([J,Re])[0]
        if np.isnan(C_T): #The interpolation returns nan if the initial points are outside range.
            C_T=0
            print('ERROR REQUESTED VALUE OUTSIDE BOUNDS FOR T')
        #Check if omega is a value or array, so we can use this function in both cases.
        if not isinstance( omega, (float, int) ):
            omega_val=omega[i]
        else:
            omega_val=omega
        T.append(C_T*atmosphere['rho']*(omega_val/(2*np.pi))**2*propeller['diameter']**4)

    return np.array(T)





if __name__ == "__main__":
    print("Run the main script")    