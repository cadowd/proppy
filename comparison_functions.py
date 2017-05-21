# -*- coding: utf-8 -*-
"""
Created on Nov 23 09:23:11 2016

Contains the functions for the relatively quick comparison of various motor
propeller combinations for a specified plane. For an explanation of roughly
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
import os
import matplotlib.pyplot as plt
import data_dealings as data_dealings
from scipy.interpolate import interp1d, interp2d, LinearNDInterpolator
#import scenarios
import numpy.polynomial.polynomial as poly
from operator import itemgetter
import itertools
from matplotlib import path, cm, rcParams, collections
import scenarios
import BLDC_model
import consumption_functions as optimised_consumption
import matplotlib.patches as mpatches
from matplotlib import _cntr as cntr
from matplotlib import _contour as contour_func
import pickle



def eff_space(U_range, motor, propeller, battery, atmosphere):
    """
    Function for calculating the efficiency and power surfaces of a given
    combination over a specified velocity range. Similar approach to that
    used in the main consumption calculation.
    """
    D=propeller['diameter']
    V_bat=battery['V']
    #Load constant models. This loads all the data we have for this combination, and is the only time we access the files.
    #propeller models
    C_Q_interp, C_T_interp, bounds=data_dealings.prop_data_funcs(propeller)
    
    Kv0=motor['Kv']*2*np.pi/60
    maxQ=100/Kv0 #Max torque at 100 amps. Thats an insane value but whatever.
    maxomega=V_bat*Kv0 #Max omega at battery volts. Quite a lot.
    Vinterp_func=None
    mot_bounds=[0,maxQ, 0, maxomega]
    I_func=None
        
    #Initial bounds to search for the rotational speed, should contain decent guesses for all possible flight RPMs
    omega_search_min=2/(bounds[3]*D)*2*np.pi #Don't bother with speeds that are highly unlikely to provide thrust, ie; outside bounds of surface
#    print(omega_search_min)
    omega_search_max=mot_bounds[3] #We also scale it by the battery voltage, as the surface max is based on a maximum of 25 volts
#    print(omega_search_max)
    
    search_res=len(U_range)
    omega_range=np.linspace(omega_search_min,omega_search_max,search_res)
    
    
    #We make a grid of rotational speeds and flight speeds representing the possible flight envelope
    X, Y = np.meshgrid(omega_range, U_range)
    #We calculate the thrust surface at every single one of these shits
    Ts = np.array([optimised_consumption.propT(x,y, C_T_interp, bounds, propeller,atmosphere) for x,y in zip(np.ravel(X), np.ravel(Y))])
    Tsurf = Ts.reshape(X.shape)
    Q_flight= np.array([optimised_consumption.propQ(x,y, C_Q_interp, bounds, propeller,atmosphere) for x,y in zip(np.ravel(X), np.ravel(Y))])
    Qsurf = Q_flight.reshape(X.shape)

    P_el=[]
    I_bat=[]
    eta_system=[]
    eta_motor=[]
    for i, omega in enumerate(X):
        V_mot, I_mot, Temp= BLDC_model.Vsurf(np.array(Qsurf[i]),np.array(omega), motor, atmosphere)
            
        #Filter out invalid results, where the motor voltage is higher than the battery voltage
        id_valid= (V_mot<=V_bat)
        P_el_val=V_mot*I_mot
        P_el_val[~id_valid]=None

        P_el.append(P_el_val)
        eta_motor.append(Qsurf[i]*omega/(P_el_val))
        eta_system.append(Tsurf[i]*Y[i]/(P_el_val))
        
        #Calculate battery current (the current used in determining motor limit)
        I_bat_val=V_mot/V_bat*I_mot
        I_bat_val[~id_valid]=None
        I_bat.append(I_bat_val)
        
#    print(Y)
#    print(P_el)
    #Make an interpolation function for repositioning the points later on
    Pelinterp_func= LinearNDInterpolator(np.transpose([np.ravel(Y), np.ravel(Tsurf)]), np.ravel(P_el), rescale=True)
    
    return Y, Tsurf, eta_system, P_el, I_bat, Pelinterp_func

def check_Psurfs(motors,propellers,Urange, Trange, folder, battery, atmosphere):
    """
    Function to check the existance of the Psurfs in a chosen folder.
    First it checks if the current settings are the same as the last time the Psurfs
    were generated. If they aren't the previous Psurfs are promptly incinerated
    and all have to be generated again. Returns a list of the merged names
    of all the Psurfs that need to be generated, as well as recreating the
    generation properties file - YOU NEED TO RUN THE PSURF GENERATION FUNCTION
    DIRECTLY AFTERWARDS.
    """
    generation_properties={
    'battery':battery,
    'atmosphere': atmosphere,
    'Urange':Urange,
    'Trange':Trange
    }
    
    regen=False
    if not os.path.exists(folder):
        os.makedirs(folder)
    try:
        last_gen_props = pickle.load( open( folder + 'generation_properties.pk', "rb" ) )
        for key in last_gen_props.keys():
            equal=last_gen_props[key]==generation_properties[key]
            if not np.all(equal):
                print('Current generation settings different to previous run. Regenerating all surfaces.')
                regen=True
    except IOError:
        pass
    
    #We err on the side of safety and instantly delete ALL the power surfaces if the settings have changed.
    #This could end up being really annoying. I'll see. There's smarter ways of doing this anyway.
    if regen:
        for root, dirs, files in os.walk(folder):
            for file in files:
                if file.endswith(".pk"):
                    os.remove(folder + file)
    combo_list=[]
    
    for motor in motors:
        for propeller in propellers:
            combo_name=motor['name'] +' - ' + propeller['name']
            save_name=folder + motor['name'] + propeller['name'] + '.pk'
            if os.path.isfile(save_name) and regen==False:
                print(combo_name + ' already calculated')
                continue
            combo_dict={
            'motor': motor,
            'propeller': propeller}
            combo_list.append(combo_dict)
        #Dump this steamer
    try:
        with open(folder + 'generation_properties.pk', 'wb') as f:
            pickle.dump(generation_properties, f)
    except Exception as e:
        print(e)
        
    return combo_list

def gen_single_Psurf(combo_dict,Urange, Trange, folder, battery, atmosphere):
    """
    Function to generate a single Psurf.
    """
    Umesh, Tmesh = np.meshgrid(Urange, Trange)
    
    motor=combo_dict['motor']
    propeller=combo_dict['propeller']
    
    combo_name=motor['name'] +' - ' + propeller['name']
    save_name=folder + motor['name'] + propeller['name'] + '.pk'

    #Here we do the actually difficult calculation. 
    print('Calculating power surface for motor: ' + motor['name'] + ' propeller: ' + propeller['name'])
    U, T, eta_system, P_el, I_bat, Pelinterp_func = eff_space(Urange, motor, propeller, battery, atmosphere)
    
    #Reinterpolate it onto the mesh. All hail the mesh.
    Psurf = Pelinterp_func(Umesh,Tmesh)
    

    #Get the maximum thrust and current in static conditions
    Imax, rpm_static, Tmax=scenarios.static_max_current(1.0,atmosphere,battery,motor,propeller)
#            Psurf_col=np.dstack((Psurf_col,Psurf)) #Collect power surface
#            motor_list.append(motor)
#            prop_list.append(propeller)
    
    results={
    'Psurf':Psurf,
    'combo': combo_name,
    'motor':motor, 
    'propeller':propeller,
    'Urange':Urange,
    'Trange':Trange,
    'Tmax':Tmax,
    'Imax': Imax,
    'Pinterp':Pelinterp_func,
    'Umesh': U,
    'Tmesh': T,
    'Pmesh': P_el,
    'eta_system': eta_system,
    }
    
    #Dump this steamer
    with open(save_name, 'wb') as f:
        pickle.dump(results, f)

def gen_Psurfs(combo_list,Urange, Trange, folder, battery, atmosphere):
    """
    Function to generate the power surface curves for each motor and propeller
    combination possible in the given lists, over the range of velocites and 
    thrusts given. The curves are saved as 2d linear interpolation functions
    into individual pickle files, which can then be re-read later as this function
    can take a long time.
    """
    Umesh, Tmesh = np.meshgrid(Urange, Trange)
#    Psurf_col=np.empty(np.shape(Umesh))
#    motor_list=[]
#    prop_list=[]
#    combo_name_list=[]
    
    generation_properties={
    'battery':battery,
    'atmosphere': atmosphere,
    'Urange':Urange,
    'Trange':Trange
    }
    
    regen=False
    try:
        last_gen_props = pickle.load( open( folder + 'generation_properties.pk', "rb" ) )
        for key in last_gen_props.keys():
            equal=last_gen_props[key]==generation_properties[key]
            if not np.all(equal):
                print('Current generation settings different to previous run. Regenerating all surfaces.')
                regen=True
    except IOError:
        pass
    #We err on the side of safety and instantly delete ALL the power surfaces if the settings have changed.
    #This could end up being really annoying. I'll see. There's smarter ways of doing this anyway.
    if regen:
        for root, dirs, files in os.walk(folder):
            for file in files:
                if file.endswith(".pk"):
                    os.remove(file)
    #Dump this steamer
    with open(folder + 'generation_properties.pk', 'wb') as f:
        pickle.dump(generation_properties, f)
    
    for combo_dict in combo_list:
        motor=combo_dict['motor']
        propeller=combo_dict['propeller']
        
        combo_name=motor['name'] +' - ' + propeller['name']
        save_name=folder + motor['name'] + propeller['name'] + '.pk'
        if os.path.isfile(save_name) and regen==False:
            print('Combination already calculated')
            continue
        print('Calculating power surface for motor: ' + motor['name'] + ' propeller: ' + propeller['name'])
        U, T, eta_system, P_el, I_bat, Pelinterp_func = eff_space(Urange, motor, propeller, battery, atmosphere)
        #Reinterpolate it onto the mesh. All hail the mesh.
        Psurf = Pelinterp_func(Umesh,Tmesh)
        

        #Get the maximum thrust and current in static conditions
        Imax, rpm_static, Tmax=scenarios.static_max_current(1.0,atmosphere,battery,motor,propeller)
        
        results={
        'Psurf':Psurf,
        'combo': combo_name,
        'motor':motor, 
        'propeller':propeller,
        'Urange':Urange,
        'Trange':Trange,
        'Tmax':Tmax,
        'Imax': Imax,
        'Pinterp':Pelinterp_func,
        'Umesh': U,
        'Tmesh': T,
        'Pmesh': P_el,
        'eta_system': eta_system,
        }
        
        #Dump this steamer
        with open(save_name, 'wb') as f:
            pickle.dump(results, f)


def calc_best_regions(motors, propellers, folder, conditions):
    """
    Function to determine the best combination for every thrust and speed
    value over a region. Once the zones in which each combination is most
    efficient is known, a small list of the best combinations for a given
    thrust speed curve (for a plane) can be VERY quickly calculated.
    """
    best_regions=[]
    combo_name_list=[]
    result_list=[]
    gen_props = pickle.load( open( folder + 'generation_properties.pk', "rb" ) )
    #Make sure the velocity and thrust ranges match, if they don't... Well it's probably not too bad but it's better to play safe.
    Umesh, Tmesh = np.meshgrid(gen_props['Urange'], gen_props['Trange'])
    Psurf_col=np.empty(np.shape(Umesh))
    #This function runs very quickly, provided all the required power surfaces
    #have already been calculated. Don't run unless
    #Check to make sure all requested power surfaces have been calculated.
    for motor in motors:
        for propeller in propellers:
            combo_name=motor['name'] +' - ' + propeller['name']
            save_name=folder + motor['name'] + propeller['name'] + '.pk'
            if os.path.isfile(save_name):
                result=pickle.load( open( save_name, "rb" ) )
                
                combo_name_list.append(combo_name)
                result_list.append(result)
                Psurf_col=np.dstack((Psurf_col,result['Psurf'])) #Collect power surface

            else:
                print('Combination missing')
                return
            
            
            
            
    #Get rid of the first entry of empties
    Psurf_col=Psurf_col[:,:,1:]
    
    for i, combo_name in enumerate(combo_name_list):
        #Now we make a merged surface of ALL of the power surfaces EXCEPT
        #the one for the current combination
        Psurf_min=np.nanmin(np.delete(Psurf_col, i, axis=2), axis=2)
#        print(np.shape(Psurfs))
        #And here we find the difference between this merged surface and the current power surface
        Diffsurf=Psurf_col[:,:,i]-Psurf_min
        #The next step is to find the 0 contour of this difference surface, this equates
        #to the intersetctions of the surfaces.

        zmin=np.nanmin(Diffsurf)
        if zmin==0:
            print(Diffsurf)
        if zmin<0:
            levels = [np.nanmin(Diffsurf),0]
            
            #Now we do some serious black magic to get the filled contour regions WITHOUT the plotting overhead
            mask=None #Don't even ask
            corner_mask=rcParams['contour.corner_mask'] #Seriously don't ask
    
            Diffsurf=np.ma.masked_invalid(Diffsurf, copy=False) #Get a masked array for the z values
            #Now we use an undocumented contour generating back-end function. May easily change in future
            #versions, watch out.
            contf= contour_func.QuadContourGenerator(Umesh, Tmesh, Diffsurf.filled(), mask, corner_mask, 0)
            #It returns a list of lists vertices and path codes (used to define the vertex type)
            vertices, kinds=contf.create_filled_contour(levels[0], levels[1])
    
            #Turn each list in the list into a path, then put the paths in another list
            paths=[path.Path(seg, codes=kind) for seg, kind in zip(vertices, kinds)]
            region={'combo_name': combo_name,
                    'region': paths}
            best_regions.append(region)
        else:
            continue
        #TODO: Get bounds of each region, and use these bounds to get an estimation
        #of possible enclosing regions for each point in order to only have to
        #check if a given point is inside 3 or 4 regions instead of all of them
        

    return best_regions

def top_choices(plane, motors, propellers, atmosphere,  folder, conditions, number):
    """
    Function that returns the best motor-prop combinations for a specified plane
    subject to a number of conditions. The conditions can be turned on or off
    and the number of top results the function returns is given by the 
    'number' parameter.
    """
    best_combos=[]
    combo_name_list=[]
    result_list=[]
    no_low_thrust=0
    no_over_current=0
    no_over_mass=0
    no_too_weak=0
    no_analysed=0
    no_suitable=0
    
    #To be safe we use the urange and trange from the generated power surface
    #properties.
    gen_props = pickle.load( open( folder + 'generation_properties.pk', "rb" ) )
    Urange=gen_props['Urange']
    Trange=gen_props['Trange']
    #Make sure the velocity and thrust ranges match, if they don't... Well it's probably not too bad but it's better to play safe.
    Umesh, Tmesh = np.meshgrid(Urange, Trange)
    Psurf_col=np.empty(np.shape(Umesh))
    
    Capacity=gen_props['battery']['capacity']*3.6*gen_props['battery']['V']
    Tlimit=conditions['static_thrust']
    masslimit=conditions['motor_mass']
    if conditions['static_thrust']==False:
        Tlimit=0
    if conditions['allowable_overload']==False:
        overload=float('inf')
    else:
        overload=conditions['allowable_overload']
    if conditions['motor_mass']==False:
        masslimit=float('inf')

    #Check to make sure all requested power surfaces have been calculated and filter them by conditions.
    for motor in motors:
        for propeller in propellers:
            no_analysed += 1
            combo_name=motor['name'] +' - ' + propeller['name']
            save_name=folder + motor['name'] + propeller['name'] + '.pk'
            if os.path.isfile(save_name):
                result=pickle.load( open( save_name, "rb" ) )
                Ilimit=overload*result['motor']['Imax']
                if result['Tmax']<Tlimit:
                    no_low_thrust +=1
                    continue
                if result['Imax']>Ilimit:
                    no_over_current +=1
                    continue
                if result['motor']['mass'] > masslimit:
                    no_over_mass +=1
                    continue
                #ADDD MAX SPEED CONDITION
                combo_name_list.append(combo_name)
                Pfunc=result['Pinterp']
                
                drag_min_vel=sp.optimize.minimize_scalar(lambda x : optimised_consumption.dragFunc(x, plane, atmosphere), bounds=(0.5, max(Urange)), method='bounded')
                d_max_vel=drag_min_vel.x
                drag_min=optimised_consumption.dragFunc(d_max_vel, plane, atmosphere)
                P_d_max=Pfunc(d_max_vel, drag_min)
                
                if P_d_max!=P_d_max:
                    no_too_weak +=1
                    continue
#                P_min_vel=sp.optimize.minimize_scalar(lambda x : Pfunc(x, optimised_consumption.dragFunc(x, plane, atmosphere)), bounds=(0.5, max(Urange)), method='bounded')
#                P_min_vel=P_min_vel.x
#                P_min=Pfunc(P_min_vel, optimised_consumption.dragFunc(P_min_vel, plane, atmosphere))
                no_suitable += 1
                if conditions['fixed_speed']==False:
                    try:
                        P_min_vel=sp.optimize.minimize_scalar(lambda x : Pfunc(x, optimised_consumption.dragFunc(x, plane, atmosphere)), bounds=(0.5, max(Urange)), method='bounded')
                        P_min_vel=P_min_vel.x
                        P_min_vel_drag=optimised_consumption.dragFunc(P_min_vel, plane, atmosphere)
                        P_min=Pfunc(P_min_vel, P_min_vel_drag)
                    except Exception as e:
                        print(e)
                        pass
                    result['max_time']=int(Capacity/P_min)
                    result['max_time_speed']=P_min_vel
                    result['max_time_drag']=P_min_vel_drag
                    result['max_time_distance']=result['max_time']*result['max_time_speed']
                    
                    result['max_distance_time']=int(Capacity/P_d_max)
                    result['max_distance_speed']=d_max_vel
                    result['max_distance_drag']=drag_min
                    result['max_distance']=result['max_distance_time']*result['max_distance_speed']
                else:
                    P_min_vel=conditions['fixed_speed']
                    P_min=Pfunc(P_min_vel, optimised_consumption.dragFunc(P_min_vel, plane, atmosphere))
                    result['max_time']=int(Capacity/P_min)
                    result['max_time_speed']=P_min_vel
                    result['max_time_distance']=result['max_time']*result['max_time_speed']

                    d_max_vel=conditions['fixed_speed']
                    P_d_max=Pfunc(d_max_vel, optimised_consumption.dragFunc(d_max_vel, plane, atmosphere))

                    result['max_distance_time']=int(Capacity/P_d_max)
                    result['max_distance_speed']=d_max_vel
                    result['max_distance']=result['max_distance_time']*result['max_distance_speed']

                
                result_list.append(result)
                Psurf_col=np.dstack((Psurf_col,result['Psurf'])) #Collect power surface

            else:
                print('Combination missing')
                return
    #Now we sort the results by whatever we're interested in
    if conditions['maximise_distance']==True:
        best_combos=sorted(result_list, key=itemgetter('max_distance'), reverse=True)[:number]
    else:
        best_combos=sorted(result_list, key=itemgetter('max_time'), reverse=True)[:number]
        
    #Summarise the analysis
    report={
    'analysed':no_analysed,
    'low_thrust':no_low_thrust,
    'over_current':no_over_current,
    'over_mass':no_over_mass,
    'too_weak':no_too_weak,
    'suitable':no_suitable,
    }
    return best_combos, report
    
    
if __name__ == "__main__":
    print("Run the main script")   
    