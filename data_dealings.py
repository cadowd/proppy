# -*- coding: utf-8 -*-
"""
Created on Nov 22 10:29:28 2016
Contains functions for dealing with the data files and extracting
approximation or interpolation functions.
The propeller data needs to be in the UI format (tests must be done in a windtunnel!)
or APC format from their website.

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
from mpl_toolkits.mplot3d import Axes3D
import csv
import scipy as sp
from itertools import islice#, izip
from scipy.interpolate import interp1d, interp2d, LinearNDInterpolator
import os
from collections import defaultdict
import numpy.polynomial.polynomial as poly
import matplotlib.cm as cm
import comparison_functions as consumption
import re




#def prop_data_static_CT(n, propeller):
#    """
#    Returns thrust coefficient of specified propeller at static conditions using UI data.
#    """
#    data_paths=[]
#    data_names=[]
#    no_files=0
#    folder=propeller['folder']
#    prop_name=propeller['name']
#    #Get names and locations of all suitable data files
#    for root, dirs, files in os.walk(folder):
#        for file in files:
#            if file.endswith(".txt") and file.startswith(prop_name):
#                data_file=os.path.join(root, file)
#                file_name=data_file.split("\\")[-1]
#                data_paths.append(data_file)
#                data_names.append(file_name)
#                no_files=no_files+1
#    
#    C_Q_points=[]
#    n_points=[]
#    C_T_points=[]
#    
#    for i,data_file in enumerate(data_paths): #
#        if 'static' in data_file:
#            with open(data_file) as csvfile:
#                reader = csv.DictReader(csvfile,skipinitialspace=True,delimiter=' ') #read the file starting from the second line
#                headers = reader.fieldnames #slurp up those nutritious header names
#                for row in reader:
#                    C_Q_points.append(float(row['CP'])/(2*np.pi))
#                    n_points.append(float(row['RPM']))
#                    C_T_points.append(float(row['CT']))
#
#                csvfile.close()
#        
#        else: 
#            pass
#            
#    
#    C_T_fit_coeffs=np.polyfit(n_points, C_T_points, 1)
#    C_T_fit = np.poly1d(C_T_fit_coeffs)
#    C_T=C_T_fit(n*60)
#    
#    return C_T

def prop_data_static_Q(omega, propeller,atmosphere):
    """
    Returns torque of specified propeller at static conditions using UI data.
    omega is in rad/s
    """
    rho=atmosphere['rho']
    D=propeller['diameter']
    if isinstance( omega, (float, int) ): #make J a list if we've fed in floats.
#        print('bloops')
        omega=np.array([omega])
    n_range=omega/(2*np.pi)
    
        
    if 'PER3' in propeller['name']: #this means APC data
        prop_data=get_APC_data(propeller)
    else: #the only other data we'll accept is UI data
        prop_data=get_UI_data(propeller)
    
    
    
    n_points=prop_data['n_points_static']
    C_Q_points=prop_data['C_Q_points_static']
    
    C_Q_interp=interp1d(n_points, C_Q_points, kind='linear')
    C_Q_fit_coeffs=np.polyfit(n_points, C_Q_points, 1)
    C_Q_fit = np.poly1d(C_Q_fit_coeffs)
    C_Q=[]
    try:
        for n in n_range:
            if n<max(n_points) and n>min(n_points):
                C_Q.append(C_Q_interp(n))
            else:
                C_Q.append(C_Q_fit(n))
    except TypeError:
        n=n_range
        if n<max(n_points) and n>min(n_points):
            C_Q.append(C_Q_interp(n))
        else:
            C_Q.append(C_Q_fit(n))
    
    if len(C_Q)>1:
        C_Q=np.array(C_Q)
    else:
        C_Q=C_Q[0]
    
    return C_Q*rho*n_range**2*D**5
    
def prop_data_static_T(omega, propeller,atmosphere):
    """
    Returns thrust of specified propeller at static conditions using UI data.
    omega is in rad/s
    """
    rho=atmosphere['rho']
    D=propeller['diameter']
    n_range=omega/(2*np.pi)
    
    if 'PER3' in propeller['name']: #this means APC data
        prop_data=get_APC_data(propeller)
    else: #the only other data we'll accept is UI data
        prop_data=get_UI_data(propeller)
    
    
    
    n_points=prop_data['n_points_static']
    C_T_points=prop_data['C_T_points_static']


    C_T_interp=interp1d(n_points, C_T_points, kind='linear')
    C_T_fit_coeffs=np.polyfit(n_points, C_T_points, 1)
    C_T_fit = np.poly1d(C_T_fit_coeffs)
    C_T=[]
    try:
        for n in n_range:
            if n<max(n_points) and n>min(n_points):
                C_T.append(C_T_interp(n))
            else:
                C_T.append(C_T_fit(n))
    except TypeError:
        n=n_range
        if n<max(n_points) and n>min(n_points):
            C_T.append(C_T_interp(n))
        else:
            C_T.append(C_T_fit(n))
    
    if len(C_T)>1:
        C_T=np.array(C_T)
    else:
        C_T=C_T[0]
#    print(n_range)
    return C_T*rho*n_range**2*D**4

def prop_data_funcs(propeller):
    """
    Returns a two dimensional linear interpolation for CQ and CT as a function of Re and J.
    Adds padding values to allow interpolation even outside of original dataset.
    """
    
    if 'PER3' in propeller['name']: #this means APC data
        prop_data=get_APC_data(propeller)
    else: #the only other data we'll accept is UI data
        prop_data=get_UI_data(propeller)
    
    
    C_Q_points=prop_data['C_Q_points']
    J_points=prop_data['J_points']
    C_T_points=prop_data['C_T_points']
    Re_points=prop_data['Re_points']
    no_points=prop_data['no_points']
    
#    n_points_static=[]
    C_T_points_static=prop_data['C_T_points_static']
    C_Q_points_static=prop_data['C_Q_points_static']
    Re_points_static=prop_data['Re_points_static']

    high_rpm_J=J_points[-no_points:] #List of the high Re values
    high_rpm_C_Q=C_Q_points[-no_points:]
    high_rpm_C_T=C_T_points[-no_points:]
    
    
    #To make out interpolation well behaved, we add a bunch of extra points that we just pull essentially out of our arses
    #This lets us build a much smoother interpolation surface and allows us to avoid any real extrapolation.
    #There is a risk this could get weirrrd for extreme values, however as long as the signs and order
    #of magnitude of our invented values are okay, the code will converge into a reasonable range, where
    #it will use values directly interpolated or ones pretty close to it.
    Re_max_data=max(Re_points)
    J_max_data=max(J_points)
    Re_max=Re_max_data*5 #Highest conceivable Re
    
    C_Q_min_data=min(C_Q_points)
    #Add static values for low J
    C_Q_points.extend(C_Q_points_static)
    C_T_points.extend(C_T_points_static)
    Re_points.extend(Re_points_static)
    J_points.extend(np.ones(len(Re_points_static))*0) #Static values start from J=0 or below
    
    #Add the static points 2 more times because the interpolation gets weird if the data isn't close together.
    C_Q_points.extend(C_Q_points_static)
    C_T_points.extend(C_T_points_static)
    Re_points.extend(Re_points_static)
    J_points.extend(np.ones(len(Re_points_static))*-0.5)
    
    #We take these values down to J=-2 (utter ridculous as a number)
    C_Q_points.extend(np.array(C_Q_points_static)*0.8) #We make it bend a little because I think it looks prettier like that
    C_T_points.extend(np.array(C_T_points_static)*0.8)
    Re_points.extend(Re_points_static)
    J_points.extend(np.ones(len(Re_points_static))*-2)
    
    #Add the equivalent 0 Re and max Re values to keep the curve nice and sexy like all the way down to 0 omega
    min_rpm_CT=C_T_points_static[0]
    min_rpm_CQ=C_Q_points_static[0]
    C_Q_points.extend([min_rpm_CQ,min_rpm_CQ,min_rpm_CQ*0.8, min_rpm_CQ,min_rpm_CQ,min_rpm_CQ*0.8])
    C_T_points.extend([min_rpm_CT,min_rpm_CT,min_rpm_CT*0.8, min_rpm_CT,min_rpm_CT,min_rpm_CT*0.8])
    Re_points.extend([0,0,0, Re_max, Re_max, Re_max])
    J_points.extend([0,-0.5,-2, 0,-0.5,-2])
    
    #Now we pad in the other direction for high J values (up to twice the max J value in the data), again a ridiculous number.
    Re_points.extend(Re_points_static) #We add our padding at around the same Re values as in the data
    J_points.extend(np.ones(len(Re_points_static))*J_max_data*1.05) #Start by adding them *just* above the highest J in the data
    C_T_points.extend(np.zeros(len(Re_points_static))) #CT is zero here, the thrust would actually be negative but we don't need that kind of stress in our lives.
    C_Q_points.extend(np.ones(len(Re_points_static))*C_Q_min_data) #CQ we force to be the minimum value.
    
    Re_points.extend(Re_points_static) #We add our padding at around the same Re values as in the data
    J_points.extend(np.ones(len(Re_points_static))*J_max_data*2.0) #Start by adding them *just* above the highest J in the data
    C_T_points.extend(np.zeros(len(Re_points_static)))
    C_Q_points.extend(np.zeros(len(Re_points_static)))
    
    #Add the equivalent 0 Re values to keep the curve nice and sexy like all the way down to 0 omega
    C_Q_points.extend([0,0, 0,0])
    C_T_points.extend([0,0, 0,0])
    Re_points.extend([0,0, Re_max, Re_max])
    J_points.extend([J_max_data*1.05,J_max_data*2.0, J_max_data*1.05,J_max_data*2.0])
    
    #Add high Re values
    Re_points.extend(np.ones(no_points)*Re_max) #We add our padding at around the same Re values as in the data
    J_points.extend(high_rpm_J) #Start by adding them *just* above the highest J in the data
    C_T_points.extend(high_rpm_C_T)
    C_Q_points.extend(high_rpm_C_Q)
    
    
    #Do a 2D linear interpolation of the surface. Can get a little weirrrrrd for far out values, but in general is okay
    CTinterp_func= LinearNDInterpolator(np.transpose([J_points, Re_points]), C_T_points, rescale=True)
    CQinterp_func= LinearNDInterpolator(np.transpose([J_points, Re_points]), C_Q_points, rescale=True)

    bounds=[0,Re_max, -2, J_max_data*2] #Surface bounds, to stop us being silly buggers
    
    return CQinterp_func, CTinterp_func, bounds#,C_T_lin_coeffs#C_Q_of_C_T_func

def get_UI_data(propeller):
    """
    Returns a two dimensional linear interpolation for CQ and CT as a function of Re and J.
    Adds padding values to allow interpolation even outside of original dataset.
    """
    
    data_paths=[]
    data_names=[]
    no_files=0
    folder=propeller['folder']
    prop_name=propeller['name']
    #Get names and locations of all suitable data files
    for root, dirs, files in os.walk(folder):
        for file in files:
            if file.endswith(".txt") and file.startswith(prop_name) and 'geom' not in file:
                data_file=os.path.join(root, file)
                file_name=data_file.split("\\")[-1]
                data_paths.append(data_file)
                data_names.append(file_name)
                no_files=no_files+1
    
    C_Q_points=[]
    J_points=[]
    C_T_points=[]
    Re_points=[]
    
    n_points_static=[]
    C_T_points_static=[]
    C_Q_points_static=[]
    Re_points_static=[]

    
    for i,data_file in enumerate(data_paths): #
        if 'static' in data_file:
            with open(data_file) as csvfile:
                reader = csv.DictReader(csvfile,skipinitialspace=True,delimiter=' ') #read the file starting from the second line
                headers = reader.fieldnames #slurp up those nutritious header names
                for row in reader:
                    C_Q_points_static.append(float(row['CP'])/(2*np.pi))
                    n_points_static.append(float(row['RPM'])/60)
                    C_T_points_static.append(float(row['CT']))
                    tip_vel_rot=float(row['RPM'])*2*np.pi/60*propeller['diameter']/2
                    Re=tip_vel_rot*propeller['diameter']/(1.5*10**-5)
                    Re_points_static.append(Re)
    
                csvfile.close()
        
        else:
            [str_start,str_end]=data_file.rsplit('_',1) #We guess that the final under score separates the name from the rpm value
            if str_end[-3:]=='txt': #If the string to the right of the final underscore is .csv, we think we've guessed right
                rpm=float(str_end[:-4])
                with open(data_file) as csvfile:
                    reader = csv.DictReader(csvfile,skipinitialspace=True,delimiter=' ') #read the file starting from the second line
                    headers = reader.fieldnames #slurp up those nutritious header names
                    no_points=0
                    for row in reader:
                        no_points=no_points+1
                        C_Q_points.append(float(row['CP'])/(2*np.pi))
                        J_points.append(float(row['J']))
                        CT=float(row['CT'])
                        #We avoid turbining data with CT<0
                        if CT<0:
                            CT=0
                        C_T_points.append(CT)
                        tip_vel_rot=rpm*2*np.pi/60*propeller['diameter']/2
                        forward_vel=float(row['J'])*rpm/60*propeller['diameter']
                        tip_vel_tot=np.sqrt(tip_vel_rot**2+forward_vel**2)
                        Re=tip_vel_tot*propeller['diameter']/(1.5*10**-5)
                        Re_points.append(Re)
                        #print(Re)
#                            plt.scatter(float(row['J']), Re, color=c)

                    csvfile.close()
    prop_data={
    'Re_points':Re_points,
    'J_points':J_points,
    'C_T_points':C_T_points,
    'C_Q_points':C_Q_points,
    'Re_points_static':Re_points_static,
    'C_T_points_static':C_T_points_static,
    'C_Q_points_static':C_Q_points_static,
    'n_points_static':n_points_static,
    'no_points':no_points}
    return prop_data

def get_APC_data(propeller):
    """
    Function to extract the data from an APC propeller data file
    """
    
    rpm_list=[]
    n_list=[]
    folder=propeller['folder']
    prop_name=propeller['name']
    data_file= folder+prop_name+".dat" #"./APC_prop_data/PER3_12x47SF.dat"

    h=[]
    C_Q_points=[]
    J_points=[]
    C_T_points=[]
    Re_points=[]
    
    C_T_points_static=[]
    C_Q_points_static=[]
    Re_points_static=[]
    n_points_static=[]

    with open(data_file) as datfile:
        '''Pass preamble'''
        n = 1 #We pass over three lines as we collect up the tasty test info goodness
        for line in datfile.readlines(): #Now we look for the rest of the data, finding the next header row
            n += 1
            if ".dat"in line:
                prop_dim_string=line.split()[0]
                D=float(prop_dim_string.rsplit('x')[0])*0.0254
            if 'PROP RPM' in line: # line with field names was found, these are the other headers
                h = [x.strip() for x in line.split('=')]
#                print(h)
                rpm_list.append(float(h[1]))
                n_list.append(n)
        
        datfile.close()
        
        for i,n in enumerate(n_list):
            csvfile = islice(open(data_file, "r"), n, None) #reopen
            reader = csv.DictReader(csvfile, skipinitialspace=True, delimiter=' ')
#            headers = reader.fieldnames
#            print(headers)
            next(reader)
            no_points=0
            rpm=rpm_list[i]
            for row in reader:
                if row['J']==None:
                    break
                if float(row['J'])==0:
                    C_Q_points_static.append(float(row['Cp'])/(2*np.pi))
                    C_T_points_static.append(float(row['Ct']))
                    tip_vel_rot=rpm*2*np.pi/60*D/2
                    Re=tip_vel_rot*D/(1.5*10**-5)
                    Re_points_static.append(Re)
                    n_points_static.append(rpm/60)
                
                no_points=no_points+1
                C_Q_points.append(float(row['Cp'])/(2*np.pi))
                J_points.append(float(row['J']))
                CT=float(row['Ct'])
                #We avoid turbining data with CT<0
                if CT<0:
                    CT=0
                C_T_points.append(CT)
                
                tip_vel_rot=rpm*2*np.pi/60*D/2
                forward_vel=float(row['J'])*rpm/60*D
                tip_vel_tot=np.sqrt(tip_vel_rot**2+forward_vel**2)
                Re=tip_vel_tot*D/(1.5*10**-5)
                Re_points.append(Re)
            
    prop_data={
    'Re_points':Re_points,
    'J_points':J_points,
    'C_T_points':C_T_points,
    'C_Q_points':C_Q_points,
    'Re_points_static':Re_points_static,
    'C_T_points_static':C_T_points_static,
    'C_Q_points_static':C_Q_points_static,
    'n_points_static':n_points_static,
    'no_points':no_points}
    return prop_data

def generate_xlfr5_constants(data_file):
    alpha_data=[]
    beta_data=[]
    CL_data=[]
    CDi_data=[]
    CDv_data=[]
    CD_data=[]
    CY_data=[]
    Cl_data=[]
    Cm_data=[]
    Cn_data=[]
    Cni_data=[]
    Qinf_data=[]
    XCP_data=[]
    h=[]
    with open(data_file) as csvfile:
        '''Pass preamble'''
        n = 3 #We pass over three lines as we collect up the tasty test info goodness
        next(csvfile)
        reader = csv.DictReader(csvfile, delimiter=',') #read the file starting from the second line
        headers = reader.fieldnames #slurp up those nutritious header names
        row=next(reader) #now we move to the next line for the filling field values
#        for (k,v) in row.items(): # go over each column name and value 
#            wing_info[k].append(v) #put each into our test_info tummy
        for line in csvfile.readlines(): #Now we look for the rest of the data, finding the next header row
            n += 1
            if 'alpha' in line: # line with field names was found, these are the other headers
                h = [x.strip() for x in line.split(',')]
                break
        csvfile.close()

        if not h==[]:
            csvfile = islice(open(data_file, "r"), n, None)
            reader = csv.DictReader(csvfile,fieldnames = h, delimiter=',')
            headers = reader.fieldnames
            try:
                for row in reader:
                    if not row['alpha']: #Break on empty row
                        break
    #                print(row)
                    alpha_data.append(float(row['alkha']))
                    beta_data.append(float(row['Beta']))
                    CL_data.append(float(row['CL']))
                    CDi_data.append(float(row['CDi']))
                    CDv_data.append(float(row['CDv']))
                    CD_data.append(float(row['CD']))
                    CY_data.append(float(row['CY']))
                    Cl_data.append(float(row['Cl']))
                    Cm_data.append(float(row['Cm']))
                    Cn_data.append(float(row['Cn']))
                    Cni_data.append(float(row['Cni']))
                    Qinf_data.append(float(row['QInf']))
                    XCP_data.append(float(row['XCP']))
            except KeyError:
                coeffs=[0,0,0]
                return coeffs
        else:
            coeffs=[0,0,0]
            return coeffs
    
#    results={
#    'y_span_data': np.array(y_span_data),
#    'chord_data': np.array(chord_data),
#    'Ai_data': np.array(Ai_data),
#    'Cl_data': np.array(Cl_data),
#    'PCd_data': np.array(PCd_data),
#    'ICd_data': np.array(ICd_data),
#    'CmGeom_data': np.array(CmGeom_data),
#    'Cm4chord_data': np.array(Cm4chord_data),
#    'XTrtop_data': np.array(XTrtop_data),
#    'XTrBot_data': np.array(XTrBot_data),
#    'XCP_data': np.array(XCP_data),
#    'BM_data': np.array(BM_data),
#    }
    coeffs=np.polyfit(CL_data, CD_data, 2)
#    print(coeffs)
#    plt.plot(CL_data,CD_data)
#    plt.show()
    return coeffs


def read_dict_file(file_name):
    """
    Reads dictionary from a file, // is a comment line, # is also a comment
    The dictionary key and value is separated on each line by a :
    It removes all trailing and leading whitespace from the values and keys,
    treats the values 'True', 'true', 'False', 'false' as booleans, floats and ints
    as floats and all else as strings.
    """
    #Regular expressions to find key and values from defaults.txt
    dict_valueRE = re.compile("(?<=:)(.*?)(?=#|\n|\Z)") #matches everything after a :
    dict_keyRE = re.compile("[^:]+") #matches everything before:
     
    component_dict = dict()
    #Reads program defaults from file into dictionary
    with open(file_name,"r") as dict_file:
            for line in dict_file:
                if not line[:2] == '//':
                    value = re.search(dict_valueRE, line).group(0).strip()
                    key = re.search(dict_keyRE, line).group(0).strip()
                    if value in ['True', 'true']:
                        value=True
                    elif value in ['False', 'false']:
                        value=False
                    else:
                        try:
                            value=float(value)
                        except ValueError:
                            pass
                    #key,value = line.strip().split(":")
                    #value.replace(" ", "")
                    component_dict[key] = value
    return component_dict
    
def create_dict_file(file_name, component_dict, header):
     
    #Reads program defaults from file into dictionary
    with open(file_name,"w") as dict_file:
        dict_file.write(header + "\n")
        for i in component_dict.keys():            
            dict_file.write(i + " : " + str(component_dict[i]) + "\n")
    return

if __name__ == "__main__":

    
    colors = iter(cm.rainbow(np.linspace(0, 1, 10)))
    data_file="./xflr5_results/T2-VLM2.csv"
    generate_xlfr5_constants(data_file)
#    get_APC_data(propeller1)