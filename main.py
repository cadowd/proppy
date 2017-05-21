# -*- coding: utf-8 -*-
"""
Interface for codes to calculate propeller and motor performance.

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

import sys
from PyQt5 import QtCore, QtGui, uic, QtWidgets

import os
from matplotlib.figure import Figure
import itertools
import numpy as np
import scipy as sp
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
    
import consumption_functions
import data_dealings
import optimise_window as optimisation
import quad_optimisation
#import estimate_plane_window as estimate_window
import pickle
import scenarios
from scipy.interpolate import interp1d
 
qtCreatorFile = "./UI/main_window.ui" # Enter file here.
 
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)


def return_checked_values(listView):
    checked_items=[]
    
    item_model = listView.model()
    for row in range(item_model.rowCount()):
        item = item_model.item(row)
        if item.checkState() == QtCore.Qt.Checked:
            checked_items.append(item.text())
    return checked_items

def return_selected_values(listView):
    checked_items=[]
    
    item_model = listView.model()
    for row in listView.selectedIndexes():
        item = item_model.itemFromIndex(row)
        checked_items.append(item.text())
    return checked_items

def get_dat_names(folder):
    #Get names and locations of all suitable data files
    file_names=[]
    for root, dirs, files in os.walk(folder):
        for file in files:
            if file.endswith(".dat"):
                data_file=os.path.join(root, file)
                file_name=data_file.split("/")[-1]
                file_name=file_name.rsplit(".",1)[0]
                file_names.append(file_name)
                
    return file_names
    
def get_prop_names(folder):
    #Get names and locations of all suitable data files
    file_names=[]
    prop_names=[]
    for root, dirs, files in os.walk(folder):
        for file in files:
            if (file.endswith(".txt") and 'static' in file):
                data_file=os.path.join(root, file)
                file_name=data_file.split("/")[-1]
                file_name=file_name.rsplit(".",1)[0]
                file_names.append(file_name)
    for name in file_names:
        [str_start,str_end]=name.rsplit('_static',1)
        prop_names.append(str_start)
    return prop_names
    
def get_prop_names_APC(folder):
    #Get names and locations of all suitable data files
    file_names=[]
    prop_names=[]
    for root, dirs, files in os.walk(folder):
        for file in files:
            if file.endswith(".dat"):
                data_file=os.path.join(root, file)
                file_name=data_file.split("/")[-1]
                file_name=file_name.rsplit(".",1)[0]
                file_names.append(file_name)
#    for name in file_names:
#        [str_start,str_end]=name.rsplit('PER3_',1)
#        prop_names.append(str_end)
    prop_names=file_names
    return prop_names

def get_prop_dict(self, prop_dat):
    #Get the props diameter from its name if it's UI data, please let it be labelled properly.
    #We also assume the diameter is given in inches as decreed by the king 400 years ago.
    if 'PER3' in prop_dat: #this is APC data boys
        folder="./APC_prop_data/"
        data_file= folder+prop_dat+".dat"
        with open(data_file) as datfile:
            n = 1
            for line in datfile.readlines(): #Now we look for the rest of the data, finding the next header row
                n += 1
                if ".dat"in line:
                    prop_dim_string=line.split()[0]
                    prop_dia=float(prop_dim_string.rsplit('x')[0])*0.0254
            datfile.close()
    else: #here be UI data
        folder="./UI_prop_data/"
        try:
            [str_start,str_end]=prop_dat.rsplit('x',1)
            prop_dia=float(str_start.rsplit('_')[-1])*0.0254
        except Exception:
            QtWidgets.QMessageBox.warning(self, 'Error in propeller file name.',
                                prop_dat + " does not match the required naming format. No diameter could be calculated.",
                                QtWidgets.QMessageBox.Ok)
            return
    propeller={
    'name':prop_dat, 
    'folder':folder, 
    'diameter': prop_dia
    }
    
    return propeller

def do_calc(self, motor_dat, prop_dat, plane_dat):
    """
    Function to do all the calculations for a given combination
    """
    battery=self.battery
    atmosphere=self.atmosphere
#    controller=self.controller

    motor=data_dealings.read_dict_file('./motors/' + motor_dat + '.dat')
    propeller=get_prop_dict(self, prop_dat)
    plane=data_dealings.read_dict_file('./planes/' + plane_dat + '.dat')
    motor_funcs, prop_funcs=consumption_functions.funky_generator(atmosphere,battery,motor,propeller)
        
    results=consumption_functions.P_el_calculator(motor_funcs,prop_funcs,plane, battery, atmosphere)
    if results==[]:
        QtWidgets.QMessageBox.warning(self, 'Not enough thrust for flight',
                    "The selection of " + propeller['name'] + " and " + motor['name'] + " is atrocious, this combination does not provide enough thrust for flight, check the battery settings perhaps.",
                    QtWidgets.QMessageBox.Ok)
    U_array=np.array([result['U'] for result in results])
    rpm_array=[result['omega']*30/np.pi for result in results]
    Pel_array=np.array([result['P_el'] for result in results])
    P_mech_prop_array=np.array([result['P_mech_prop'] for result in results])
    P_mech_shaft_array=np.array([result['P_mech_shaft'] for result in results])
    
    Th_array=P_mech_prop_array/U_array
    eta_mot=P_mech_shaft_array/Pel_array
    eta_prop=P_mech_prop_array/P_mech_shaft_array
    eta_total=eta_mot*eta_prop
    
    line_label=plane_dat + ' ' + motor_dat + ' ' + prop_dat

    #Find the location of minimum drag
    drag_min_vel=sp.optimize.minimize_scalar(lambda x : consumption_functions.dragFunc(x, plane, atmosphere), bounds=(0.5, 100), method='bounded')
    drag_min_vel=drag_min_vel.x
    
#    drag_min=consumption_functions.dragFunc(drag_min_vel, plane, atmosphere)

    result={
    'line_label': line_label,
    'U': U_array,
    'Pel': Pel_array,
    'rpm': rpm_array,
    'P_mech_prop': P_mech_prop_array,
    'P_mech_shaft': P_mech_shaft_array,
    'Th': Th_array,
    'eta_mot': eta_mot,
    'eta_prop': eta_prop,
    'eta_total': eta_total,
    'drag_min_vel': drag_min_vel
    }
    return result
    
    
class main_window(QtWidgets.QMainWindow, Ui_MainWindow):
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)

        self.setupUi(self)
        
        fig1 = Figure()
        ax1f1 = fig1.add_subplot(111)
        self.canvas = FigureCanvas(fig1)
        self.layout_plot.addWidget(self.canvas)
        self.toolbar = NavigationToolbar(self.canvas, 
                self, coordinates=True)
#        self.addToolBar(self.toolbar)
#        self.axes=ax1f1
        self.layout_plot.addWidget(self.toolbar)
        self.plot_colors = itertools.cycle(['k', "r", "b", "g", "y", "c","m"])
        self.current_plots = []
        self.current_results = []
        
        #Get list of allll the propellers for which we have data
        prop_list = QtGui.QStandardItemModel(self.listView_props)
        self.full_prop_name_list=get_prop_names('./UI_prop_data/')+get_prop_names_APC('./APC_prop_data/')
        for name in self.full_prop_name_list:
            item = QtGui.QStandardItem(name)
            item.setCheckable(True)
            prop_list.appendRow(item)
        self.listView_props.setModel(prop_list)
        
        self.pop_list_views()
        
        #Menu actions        
        self.actionQuit.triggered.connect(self.close)
        self.actionBattery_settings.triggered.connect(self.battery_settings)
        self.actionAtmosphere_settings.triggered.connect(self.atmosphere_settings)
#        self.actionController_settings.triggered.connect(self.controller_settings)
        self.actionOptimisation_settings.triggered.connect(self.optimisation_settings)
                
        #Button functions and so on
        self.pushButton_plot.clicked.connect(lambda: self.plot(ax1f1))
        self.pushButton_clear_plot.clicked.connect(lambda: self.clear_plot(ax1f1))
        
        self.radioButton_V_Pel.toggled.connect(lambda: self.change_plot(ax1f1))
        self.radioButton_V_eta.toggled.connect(lambda: self.change_plot(ax1f1))
        self.radioButton_V_eta_prop.toggled.connect(lambda: self.change_plot(ax1f1))
        self.radioButton_V_eta_mot.toggled.connect(lambda: self.change_plot(ax1f1))
        
        
        self.pushButton_motor_edit.clicked.connect(self.edit_motor)
        self.pushButton_motor_create.clicked.connect(self.new_motor)
        
        self.pushButton_plane_edit.clicked.connect(self.edit_plane)
        self.pushButton_plane_create.clicked.connect(self.create_plane)
        
        self.pushButton_prop_filter.clicked.connect(self.filter_props)
        
        
        self.pushButton_optimise.clicked.connect(self.optimise_window)
        self.pushButton_quad_optimise.clicked.connect(self.quad_optimise_window)
        self.pushButton_get_risks.clicked.connect(self.get_risks)
        self.pushButton_prop_dets.clicked.connect(self.kv_functions)
        
        self.opt_settings=pickle.load( open( 'opt_settings.pk', "rb" ) )

        
        self.battery={
        'V':14.8, #battery voltage
        'capacity':4800
        }
        
        self.atmosphere={
        'rho':1.20, #
        'mu':1.8e-5, #dynamic viscosity
        'Ta': 25, #Ambient temperature
        }
        
        
    def pop_list_views(self):
        """
        Function to populate lists with the files found in the folders.
        """
        motor_list = QtGui.QStandardItemModel(self.listView_motors)
        for name in get_dat_names('./motors/'):
            item = QtGui.QStandardItem(name)
            item.setCheckable(True)
            motor_list.appendRow(item)
        self.listView_motors.setModel(motor_list)
        
        plane_list = QtGui.QStandardItemModel(self.listView_planes)
        for name in get_dat_names('./planes/'):
            item = QtGui.QStandardItem(name)
            item.setCheckable(True)
            plane_list.appendRow(item)
        self.listView_planes.setModel(plane_list)
        
#    def addmpl(self, fig):
    
    def filter_props(self):
        """
        Function to filter props by name.
        """
        filtered_prop_name_list=[]
        filter_string=self.lineEdit_prop_filter.text()
        for prop_name in self.full_prop_name_list:
            if filter_string in prop_name:
                filtered_prop_name_list.append(prop_name)
                
        prop_list = QtGui.QStandardItemModel(self.listView_props)
        for name in filtered_prop_name_list:
            item = QtGui.QStandardItem(name)
            item.setCheckable(True)
            prop_list.appendRow(item)
        self.listView_props.setModel(prop_list)
        
    def get_risks(self):
        """
        Function to calculate properties of selected (NOT CHECKED) combination.
        """
        try:
            motor_dat=return_selected_values(self.listView_motors)[0]
            prop_dat=return_selected_values(self.listView_props)[0]
            plane_dat=return_selected_values(self.listView_planes)[0]
        except Exception:
            QtWidgets.QMessageBox.warning(self, 'Select a motor, propeller and plane.',
                                            "The system of interest must be selected (different to being checked).",
                                            QtWidgets.QMessageBox.Ok)
            return        
        
        result = do_calc(self, motor_dat, prop_dat, plane_dat)
        umax=max(result['U'])
        umin=min(result['U'])
        etamax=max(result['eta_total'])*100
        I_flight=result['Pel']/self.battery['V']
        Imax_flight=max(I_flight)
        cruise_vel=result['drag_min_vel']
        motor=data_dealings.read_dict_file('./motors/' + motor_dat + '.dat')
        propeller=get_prop_dict(self, prop_dat)
        Imax_static, rpm_static, T_static, Q_static =scenarios.static_max_current(1.0,self.atmosphere,self.battery,motor,propeller)
        
        Imax_flight_percent = Imax_flight/motor['Imax']*100
        Imax_static_percent = Imax_static/motor['Imax']*100
        
        Iflight_str='{:0.2f} A ({:0.0f}% of motor limit.)'.format(Imax_flight,Imax_flight_percent)
        Istatic_str='{:0.2f} A ({:0.0f}% of motor limit.)'.format(Imax_static,Imax_static_percent)
        
        plane=data_dealings.read_dict_file('./planes/' + plane_dat + '.dat')
        
        max_acc=T_static/plane['mass']
        risk_string="""Motor: {} \n
Propeller: {} \n
Plane: {} \n
Maximum system efficiency: {:0.2f}%\n
Maximum flight speed: {:0.2f} m/s \n
Cruising flight speed: {:0.2f} m/s \n
Minimum flight speed: {:0.2f} m/s *Plane may stall above this speed\n
Maximum current in flight: {} \n
Maximum current at takeoff: {} \n
Maximum thrust at takeoff: {:0.2f} N \n
Maximum acceleration at takeoff: {:0.2f} m/s/s \n
""".format(motor_dat,prop_dat,plane_dat,etamax,umax, cruise_vel, umin,Iflight_str, Istatic_str,T_static, max_acc)

        QtWidgets.QMessageBox.information(self, 'Details of selection.',
                                risk_string,
                                QtWidgets.QMessageBox.Ok)
        
    
    def plot(self, axes):
        """
        Function to plot all selected combinations. Please don't be silly and
        try to plot everything, there's better ways to waste your time.
        """
        motors=return_checked_values(self.listView_motors)
        props=return_checked_values(self.listView_props)
        planes=return_checked_values(self.listView_planes)
        
        if len(motors)==0 or len(props)==0 or len(planes)==0:
            QtWidgets.QMessageBox.warning(self, 'Select a motor, propeller and plane.',
                                            "A motor, propeller and plane must be selected to plot.",
                                            QtWidgets.QMessageBox.Ok)
            return
        
        no_to_plot=0
        calc_list=[]
        for motor_dat in motors:
            for prop_dat in props:
                for plane_dat in planes:
                    line_label=plane_dat + ' ' + motor_dat + ' ' + prop_dat
                    if line_label in self.current_plots:
                        continue
                    no_to_plot +=1
                    calc_dict={
                    'motor':motor_dat,
                    'prop':prop_dat,
                    'plane':plane_dat,}
                    calc_list.append(calc_dict)

        if len(calc_list)==0:
            return
        self.prog_widget = QtWidgets.QProgressDialog("Plotting selection", "Stop the madness!",0,no_to_plot)
        self.prog_widget.setWindowModality(QtCore.Qt.WindowModal)
#        self.prog_widget.setMinimumDuration(0)
        self.prog_widget.show()


        # Setting the value on every run to 0
        self.prog_widget.setValue(0)

        for calc_dict in calc_list:
            QtWidgets.QApplication.instance().processEvents()
            if self.prog_widget.wasCanceled():
                QtWidgets.QMessageBox.information(self, 'Plotting cancelled',
                                        "The plotting has been cancelled.",
                                        QtWidgets.QMessageBox.Ok)
                return

            result = do_calc(self, calc_dict['motor'], calc_dict['prop'], calc_dict['plane'])
            self.prog_widget.setValue(self.prog_widget.value()+1)
            self.single_plot(result, axes)



    def single_plot(self, result, axes):
        """
        Function to draw a single plot
        """
#        print('Trying')
        line_label=result['line_label']
        if line_label in self.current_plots:
            return
        c=next(self.plot_colors)
        
        self.current_plots.append(line_label)
        self.current_results.append(result)
        
        if self.radioButton_V_Pel.isChecked():
            axes.plot(result['U'], result['Pel'], label=line_label, color=c, marker='x')
        elif self.radioButton_V_eta.isChecked():
            axes.plot(result['U'], result['eta_total'], label=line_label, color=c, marker='x')
        elif self.radioButton_V_eta_prop.isChecked():
            axes.plot(result['U'], result['eta_prop'], label=line_label, color=c, marker='x')
        elif self.radioButton_V_eta_mot.isChecked():
            axes.plot(result['U'], result['eta_mot'], label=line_label, color=c, marker='x')
#        self.addmpl(fig1)
        axes.legend(loc=1)
#        axes.axis([11, 24, 0, 200])
        if self.radioButton_V_Pel.isChecked():
            axes.set_ylabel(r"Consumed power [W]")
            axes.set_title("Power vs air speed")
        elif self.radioButton_V_eta.isChecked():
            axes.set_ylabel(r"Total system efficiency")
            axes.set_title("Efficiency vs air speed")
        elif self.radioButton_V_eta_prop.isChecked():
            axes.set_ylabel(r"Propeller efficiency")
            axes.set_title("Efficiency vs air speed")
        elif self.radioButton_V_eta_mot.isChecked():
            axes.set_ylabel(r"Motor efficiency")
            axes.set_title("Efficiency vs air speed")
        
        axes.set_xlabel(r"Air speed (in steady flight) [m/s]")
        self.canvas.draw()
    
    def clear_plot(self, axes):
        """
        Function to clear plot, resetting colours and current plot data
        """
        axes.clear()
        self.current_plots=[]
        self.plot_colors = itertools.cycle(['k', "r", "b", "g", "y", "c","m"])
        self.canvas.draw()
        
    def change_plot(self, axes):
        """
        Function to redraw plot
        """
        axes.clear()
        self.plot_colors = itertools.cycle(['k', "r", "b", "g", "y", "c","m"])
        
        for i, line_label in enumerate(self.current_plots):
            c=next(self.plot_colors)
            result=self.current_results[i]
            if self.radioButton_V_Pel.isChecked():
                axes.plot(result['U'], result['Pel'], label=line_label, color=c, marker='x')
            elif self.radioButton_V_eta.isChecked():
                axes.plot(result['U'], result['eta_total'], label=line_label, color=c, marker='x')
            elif self.radioButton_V_eta_prop.isChecked():
                axes.plot(result['U'], result['eta_prop'], label=line_label, color=c, marker='x')
            elif self.radioButton_V_eta_mot.isChecked():
                axes.plot(result['U'], result['eta_mot'], label=line_label, color=c, marker='x')
#        self.addmpl(fig1)
        axes.legend(loc=4)
#        axes.axis([11, 24, 0, 200])
        if self.radioButton_V_Pel.isChecked():
            axes.set_ylabel(r"Consumed power [W]")
            axes.set_title("Power vs air speed")
        elif self.radioButton_V_eta.isChecked():
            axes.set_ylabel(r"Total system efficiency")
            axes.set_title("Efficiency vs air speed")
        elif self.radioButton_V_eta_prop.isChecked():
            axes.set_ylabel(r"Propeller efficiency")
            axes.set_title("Efficiency vs air speed")
        elif self.radioButton_V_eta_mot.isChecked():
            axes.set_ylabel(r"Motor efficiency")
            axes.set_title("Efficiency vs air speed")
        
        axes.set_xlabel(r"Air speed (in steady flight) [m/s]")
        self.canvas.draw()
        
    def edit_motor(self):
        self.edit_motor_dialog = edit_motor_window()

        motor_dat=return_selected_values(self.listView_motors)
        if len(motor_dat)==0:
            QtWidgets.QMessageBox.warning(self, 'No motor selected to edit.',
                                            "Please select a motor from the list for editing.",
                                            QtWidgets.QMessageBox.Ok)
            return
        motor=data_dealings.read_dict_file('./motors/' + motor_dat[0] + '.dat')
        self.edit_motor_dialog.motor_edited=motor
        edit_motor_window.populate(self.edit_motor_dialog, motor, motor_dat[0])
        
        
        edit_motor_window.refresh(self.edit_motor_dialog)
        
        # For Modal dialogs
        self.edit_motor_dialog.exec_()
        self.pop_list_views()
        
        # Or for modeless dialogs
        # self.dialog.show()
#        return
        
    def new_motor(self):
        self.new_motor_dialog = new_motor_window()
        
        edit_motor_window.refresh(self.new_motor_dialog)
        
        # For Modal dialogs
        self.new_motor_dialog.exec_()
        self.pop_list_views()
        
    def edit_plane(self):
        self.plane_dialog = edit_plane_window()
        
        plane_dat=return_selected_values(self.listView_planes)
        if len(plane_dat)==0:
            QtWidgets.QMessageBox.warning(self, 'No plane selected to edit.',
                                            "Please select a plane from the list for editing.",
                                            QtWidgets.QMessageBox.Ok)
            return
        plane=data_dealings.read_dict_file('./planes/' + plane_dat[0] + '.dat')
        self.plane_dialog.plane_edited=plane_dat[0]
        edit_plane_window.populate(self.plane_dialog, plane, plane_dat[0])
        
        # For Modal dialogs
        self.plane_dialog.exec_()
        self.pop_list_views()
    
    def create_plane(self):
        self.plane_dialog = create_plane_window()
        
#        self.plane_dialog.pushButton_estimate.clicked.connect(self.estimate_plane)

        # For Modal dialogs
        self.plane_dialog.exec_()
        self.pop_list_views()
        
#    def estimate_plane(self):
#        self.dialog = estimate_window.estimate_plane_window(self.atmosphere)
#        self.dialog.pushButton_accept.clicked.connect(lambda: self.estimate_return(self.dialog))
#        
#        # For Modal dialogs
#        self.dialog.exec_()
        
    def estimate_return(self, window):
        try:
            coeffs=np.polyfit(window.Cl, window.Cd, 2)
        except AttributeError:
            QtWidgets.QMessageBox.warning(self, 'Polar has not yet been calculated',
                                            "Please fill in all boxes to calculate the polar coefficients.",
                                            QtWidgets.QMessageBox.Ok)
            return
        
        self.plane_dialog.doubleSpinBox_C1.setValue(coeffs[0])
        self.plane_dialog.doubleSpinBox_C2.setValue(coeffs[1])
        self.plane_dialog.doubleSpinBox_C3.setValue(coeffs[2])
        
        self.plane_dialog.doubleSpinBox_ref_area.setValue(window.S_ref)
        self.plane_dialog.doubleSpinBox_mass.setValue(window.total_mass)
        window.close()

    def optimisation_settings(self):
        self.dialog = optimisation_settings_window()
        self.dialog.buttonBox.accepted.connect(lambda: self.optimisation_settings_okay(self.dialog))
        optimisation_settings_window.populate(self.dialog, self.opt_settings)
        # For Modal dialogs
        self.dialog.exec_()
    
    def optimisation_settings_okay(self, dialog):
        self.opt_settings={
        'Umin':dialog.doubleSpinBox_Umin.value(),
        'Umax':dialog.doubleSpinBox_Umax.value(),
        'Tmin': dialog.doubleSpinBox_Tmin.value(),
        'Tmax': dialog.doubleSpinBox_Tmax.value(),
        'resolution':dialog.spinBox_res.value()
        }
        with open('opt_settings.pk', 'wb') as f:
            pickle.dump(self.opt_settings, f)

    def battery_settings(self):
        self.dialog = battery_settings_window()
        self.dialog.buttonBox.accepted.connect(lambda: self.battery_settings_okay(self.dialog))
        battery_settings_window.populate(self.dialog, self.battery)
        # For Modal dialogs
        self.dialog.exec_()
    
    def battery_settings_okay(self, dialog):
        self.battery={
        'V': dialog.doubleSpinBox_voltage.value(),
        'capacity': dialog.doubleSpinBox_capacity.value(),
        }
        
    def atmosphere_settings(self):
        self.dialog = atmosphere_settings_window()
        self.dialog.buttonBox.accepted.connect(lambda: self.atmosphere_settings_okay(self.dialog))
        atmosphere_settings_window.populate(self.dialog, self.atmosphere)
        # For Modal dialogs
        self.dialog.exec_()
    
    def atmosphere_settings_okay(self, dialog):
        self.atmosphere={
        'rho': dialog.doubleSpinBox_rho.value(),
        'mu': dialog.doubleSpinBox_mu.value()/10**5,
        'Ta': dialog.doubleSpinBox_Ta.value(),
        }
        
        
    def optimise_window(self):
        motors=return_checked_values(self.listView_motors)
        propellers=return_checked_values(self.listView_props)
        planes=return_checked_values(self.listView_planes)

        if len(motors)==0 or len(propellers)==0 or len(planes)==0:
            QtWidgets.QMessageBox.warning(self, 'Select a motor, propeller and plane.',
                                            "A motor, propeller and plane must be selected for optimisation.",
                                            QtWidgets.QMessageBox.Ok)
            return
            
        if len(planes)>1:
            QtWidgets.QMessageBox.warning(self, 'Too many planes selected',
                                            "Only a single plane may be selected for optimisation.",
                                            QtWidgets.QMessageBox.Ok)
            return
        
        opt_settings=self.opt_settings
        
        options={
        'battery': self.battery,
        'atmosphere': self.atmosphere,
        'optimisation': opt_settings,
        }
        self.window = optimisation.optimise_window(planes[0], propellers, motors, options)
        self.window.pushButton_return.clicked.connect(lambda: self.optimise_return(self.window))
        # For Modal dialogs
        self.window.exec_()
        
        
        
    def optimise_return(self, window):

        results_dat=[]
        result_motors=[]
        result_props=[]
    
        item_model = window.listView_results.model()
        for row in range(item_model.rowCount()):
            item = item_model.item(row)
            if item.checkState() == QtCore.Qt.Checked:
                results_dat.append(item.combo)
                result_motors.append(item.combo['motor']['name'])
                result_props.append(item.combo['propeller']['name'])
        

        
        item_model = self.listView_motors.model()
        for row in range(item_model.rowCount()):
            item = item_model.item(row)
            item.setCheckState(QtCore.Qt.Unchecked)
            if item.text() in result_motors:
                item.setCheckState(QtCore.Qt.Checked)
        self.listView_motors.setModel(item_model)
        
        item_model = self.listView_props.model()
        for row in range(item_model.rowCount()):
            item = item_model.item(row)
            item.setCheckState(QtCore.Qt.Unchecked)
            if item.text() in result_props:
                item.setCheckState(QtCore.Qt.Checked)
        self.listView_props.setModel(item_model)
        
        window.close()
        
    def quad_optimise_window(self):
        motors=return_checked_values(self.listView_motors)
        propellers=return_checked_values(self.listView_props)
#        planes=return_checked_values(self.listView_planes)

        if len(motors)==0 or len(propellers)==0:
            QtWidgets.QMessageBox.warning(self, 'Select at least one motor and one propeller',
                                            "A motor and propeller must be selected for hover craft optimisation.",
                                            QtWidgets.QMessageBox.Ok)
            return
        
        opt_settings=self.opt_settings
        
        options={
        'battery': self.battery,
        'atmosphere': self.atmosphere,
        'optimisation': opt_settings,
        }
        self.window = quad_optimisation.quad_optimise_window(propellers, motors, options)
#        self.window.pushButton_return.clicked.connect(lambda: self.optimise_return(self.window))
        # For Modal dialogs
        self.window.exec_()
        
    def propeller_functions(self):
        
        try:
            prop_dat=return_selected_values(self.listView_props)[0]
        except IndexError:
            QtWidgets.QMessageBox.warning(self, 'No propeller selected.',
                                            "Please select a propeller from the list for analysis.",
                                            QtWidgets.QMessageBox.Ok)
            return
        
        propeller=get_prop_dict(self, prop_dat)
        
        C_Q_interp, C_T_interp, bounds=data_dealings.prop_data_funcs(propeller)
        prop_funcs={'C_Q_interp':C_Q_interp,
                     'C_T_interp':C_T_interp,
                     'bounds':bounds,
                     }
        
#        self.propeller_dialog.propeller=propeller
        #get functions
        self.propeller_dialog = propeller_dialog_window(propeller,prop_funcs, self.atmosphere, self.battery)
        propeller_dialog_window.populate(self.propeller_dialog, propeller)
        
        # For Modal dialogs
        self.propeller_dialog.exec_()
        
    def kv_functions(self):
        
#        try:
#            prop_dat=return_selected_values(self.listView_props)[0]
#        except IndexError:
#            QtWidgets.QMessageBox.warning(self, 'No propeller selected.',
#                                            "Please select a propeller from the list for analysis.",
#                                            QtWidgets.QMessageBox.Ok)
#            return
        
#        propeller=get_prop_dict(self, prop_dat)
#        
#        C_Q_interp, C_T_interp, bounds=data_dealings.prop_data_funcs(propeller)
#        prop_funcs={'C_Q_interp':C_Q_interp,
#                     'C_T_interp':C_T_interp,
#                     'bounds':bounds,
#                     }
        
#        self.propeller_dialog.propeller=propeller
        #get functions
        self.kv_dialog = kv_dialog_window(self.atmosphere, self.battery)
        kv_dialog_window.populate(self.kv_dialog)
        
        # For Modal dialogs
        self.kv_dialog.exec_()


class atmosphere_settings_window(QtWidgets.QDialog):
    def __init__(self):
        super(atmosphere_settings_window, self).__init__()
        uic.loadUi("./UI/dialog_atmosphere_settings.ui", self)
        self.pushButton_calc_props.clicked.connect(self.calc_props)
    
    def populate(self, atmosphere):
        self.doubleSpinBox_rho.setValue(atmosphere['rho'])
        self.doubleSpinBox_mu.setValue(atmosphere['mu']*10**5)
        self.doubleSpinBox_Ta.setValue(atmosphere['Ta'])
        
    def calc_props(self, atmosphere):
        R_d=287.058 #Specific gas constant for dry air
        R_vap=461.495 #Specific gas constant for water vapour
        
        
        Temp=self.doubleSpinBox_Ta.value()
        humid=self.doubleSpinBox_humidity.value()
        h=self.doubleSpinBox_alt.value()
        
        p_sat=6.1708*100*10**(Temp*7.5/(Temp+273.3)) #saturation pressure
        
        p_vap=humid/100*p_sat #vapour pressure
        
        p_abs=101325*(1-2.25577e-5*h)**5.25588 
        
        
        p_d=p_abs-p_vap
        rho=p_d/(R_d*(Temp+273.3))+p_vap/(R_vap*(Temp+273.3))        
        
        self.doubleSpinBox_rho.setValue(rho)

class optimisation_settings_window(QtWidgets.QDialog):
    def __init__(self):
        super(optimisation_settings_window, self).__init__()
        uic.loadUi("./UI/dialog_optimisation_settings.ui", self)
        self.pushButton_reset_Psurfs.clicked.connect(lambda: optimisation_settings_window.reset_Psurfs(self))
    
    def populate(self, opt_settings):
        self.doubleSpinBox_Umin.setValue(opt_settings['Umin'])
        self.doubleSpinBox_Umax.setValue(opt_settings['Umax'])
        self.doubleSpinBox_Tmin.setValue(opt_settings['Tmin'])
        self.doubleSpinBox_Tmax.setValue(opt_settings['Tmax'])
        self.spinBox_res.setValue(opt_settings['resolution'])
        
    def reset_Psurfs(self):
        folder='./power_surface_precalcs/'
        for root, dirs, files in os.walk(folder):
            for file in files:
                if file.endswith(".pk"):
                    os.remove(folder + file)
        QtWidgets.QMessageBox.information(self, 'Power surfaces deleted',
                                            "All previously generated power surfaces have been deleted.",
                                            QtWidgets.QMessageBox.Ok)
        
class battery_settings_window(QtWidgets.QDialog):
    def __init__(self):
        super(battery_settings_window, self).__init__()
        uic.loadUi("./UI/dialog_battery_settings.ui", self)
    
    def populate(self, battery):
        self.doubleSpinBox_voltage.setValue(battery['V'])
        self.doubleSpinBox_capacity.setValue(battery['capacity'])

class create_plane_window(QtWidgets.QDialog):
    def __init__(self):
        super(create_plane_window, self).__init__()
        uic.loadUi("./UI/dialog_create_plane.ui", self)
        self.pushButton_xflr5.clicked.connect(lambda: edit_plane_window.get_xflr5_file(self))
        self.buttonBox.accepted.connect(lambda: edit_plane_window.okay(self))
      
class edit_plane_window(QtWidgets.QDialog):
    def __init__(self):
        super(edit_plane_window, self).__init__()
        uic.loadUi("./UI/dialog_edit_plane.ui", self)
        
        self.pushButton_xflr5.clicked.connect(lambda: edit_plane_window.get_xflr5_file(self))
        self.pushButton_new.clicked.connect(lambda: edit_plane_window.new_plane(self))
        self.buttonBox.accepted.connect(lambda: edit_plane_window.okay(self))

    def populate(self, plane, name):
        self.lineEdit_display_name.setText(name)
        self.doubleSpinBox_ref_area.setValue(plane['ref_area'])
        self.doubleSpinBox_mass.setValue(plane['mass'])
        self.doubleSpinBox_C1.setValue(plane['C_1'])
        self.doubleSpinBox_C2.setValue(plane['C_2'])
        self.doubleSpinBox_C3.setValue(plane['C_3'])
        
    def okay(self):
        file_name=self.lineEdit_display_name.text()
        plane={
        'ref_area': self.doubleSpinBox_ref_area.value(),
        'mass': self.doubleSpinBox_mass.value(),
        'C_1': self.doubleSpinBox_C1.value(),
        'C_2': self.doubleSpinBox_C2.value(),
        'C_3': self.doubleSpinBox_C3.value(),
        }
        
        try:   
            original_file_name=self.plane_edited
            try:
                os.remove('./planes/' + original_file_name + '.dat')
            except OSError:
                pass
        except AttributeError:
            pass
        
        data_dealings.create_dict_file('./planes/' + file_name + '.dat', plane, "// File defining plane constants and values to be read into dictionary")
        main_window.pop_list_views(main_window())
        
    def new_plane(self):
        file_name=self.lineEdit_display_name.text()
        motor={
        'ref_area': self.doubleSpinBox_ref_area.value(),
        'mass': self.doubleSpinBox_mass.value(),
        'C_1': self.doubleSpinBox_C1.value(),
        'C_2': self.doubleSpinBox_C2.value(),
        'C_3': self.doubleSpinBox_C3.value(),
        }
        
        data_dealings.create_dict_file('./planes/' + file_name + '.dat', motor, "// File defining plane constants and values to be read into dictionary")
        main_window.pop_list_views(main_window())
        self.close()
    def get_xflr5_file(self):
        cwd = os.getcwd()
        folder=cwd+"/xflr5_results/"
        xflr5_file = QtWidgets.QFileDialog.getOpenFileName(self, 'Open xflr5 results polar file', 
         folder,"XFLR5 polars (*.csv)")
         
        coeffs=data_dealings.generate_xlfr5_constants(xflr5_file)
        if all([ v == 0 for v in coeffs ]):
            QtWidgets.QMessageBox.warning(self, 'No coefficients extracted',
                                            "Please check to ensure file selected is a valid XFLR5 type 2 exported polar.",
                                            QtWidgets.QMessageBox.Ok)
        self.doubleSpinBox_C1.setValue(coeffs[0])
        self.doubleSpinBox_C2.setValue(coeffs[1])
        self.doubleSpinBox_C3.setValue(coeffs[2])

        
class new_motor_window(QtWidgets.QDialog):
    def __init__(self):
        super(new_motor_window, self).__init__()
        uic.loadUi("./UI/dialog_new_motor.ui", self)
                
        self.checkBox_use_thermal.clicked.connect(lambda: edit_motor_window.refresh(self))
        self.pushButton_accept.clicked.connect(lambda: edit_motor_window.okay(self))
        self.pushButton_reject.clicked.connect(lambda: edit_motor_window.cancel(self))
        
class edit_motor_window(QtWidgets.QDialog):
    def __init__(self):
        super(edit_motor_window, self).__init__()
        uic.loadUi("./UI/dialog_edit_motor.ui", self)
        
        self.pushButton_new.clicked.connect(lambda: edit_motor_window.new_motor(self))
        
        self.checkBox_use_thermal.clicked.connect(lambda: edit_motor_window.refresh(self))
        self.pushButton_accept.clicked.connect(lambda: edit_motor_window.okay(self))
        self.pushButton_reject.clicked.connect(lambda: edit_motor_window.cancel(self))
        
    def populate(self, motor, name):
        self.lineEdit_display_name.setText(name)
        self.doubleSpinBox_mass.setValue(motor['mass'])
        self.checkBox_use_thermal.setChecked(motor['use_thermal'])
        self.doubleSpinBox_Kv.setValue(motor['Kv'])
        self.doubleSpinBox_I0.setValue(motor['I0'])
        self.doubleSpinBox_Imax.setValue(motor['Imax'])
        self.doubleSpinBox_R.setValue(motor['R'])
        self.doubleSpinBox_Tref.setValue(motor['Tref'])
        self.doubleSpinBox_k1.setValue(motor['k1'])
        self.doubleSpinBox_k2.setValue(motor['k2'])
        self.doubleSpinBox_k3.setValue(motor['k3'])
        self.doubleSpinBox_Rthwm.setValue(motor['Rthwm'])
        self.doubleSpinBox_Rthma.setValue(motor['Rthma'])
        self.doubleSpinBox_Rthws.setValue(motor['Rthws'])
        self.doubleSpinBox_Rthsa.setValue(motor['Rthsa'])
        self.doubleSpinBox_alphar.setValue(motor['alphar'])
        self.doubleSpinBox_alpham.setValue(motor['alpham'])
    
    def refresh(self):            
        if (self.checkBox_use_thermal.isChecked()):
            self.doubleSpinBox_Rthwm.setEnabled(True)
            self.doubleSpinBox_Rthma.setEnabled(True)
            self.doubleSpinBox_Rthws.setEnabled(True)
            self.doubleSpinBox_Rthsa.setEnabled(True)
            self.doubleSpinBox_alphar.setEnabled(True)
            self.doubleSpinBox_alpham.setEnabled(True)
        else:
            self.doubleSpinBox_Rthwm.setEnabled(False)
            self.doubleSpinBox_Rthma.setEnabled(False)
            self.doubleSpinBox_Rthws.setEnabled(False)
            self.doubleSpinBox_Rthsa.setEnabled(False)
            self.doubleSpinBox_alphar.setEnabled(False)
            self.doubleSpinBox_alpham.setEnabled(False)
    
    def new_motor(self):
        file_name=self.lineEdit_display_name.text()
        motor={
        'use_thermal': self.checkBox_use_thermal.isChecked(),
        'name': self.lineEdit_display_name.text(),
        'mass': self.doubleSpinBox_mass.value(),
        'Kv': self.doubleSpinBox_Kv.value(),
        'I0': self.doubleSpinBox_I0.value(),
        'Imax': self.doubleSpinBox_Imax.value(),
        'R': self.doubleSpinBox_R.value(),
        'Tref': self.doubleSpinBox_Tref.value(),
        'k1': self.doubleSpinBox_k1.value(),
        'k2': self.doubleSpinBox_k2.value(),
        'k3': self.doubleSpinBox_k3.value(),
        'Rthwm': self.doubleSpinBox_Rthwm.value(),
        'Rthma': self.doubleSpinBox_Rthma.value(),
        'Rthws': self.doubleSpinBox_Rthws.value(),
        'Rthsa': self.doubleSpinBox_Rthsa.value(),
        'alphar': self.doubleSpinBox_alphar.value(),
        'alpham': self.doubleSpinBox_alpham.value(),
        }
        
        
        
        data_dealings.create_dict_file('./motors/' + file_name + '.dat', motor, "// File defining motor constants and values to be read into dictionary")
        main_window.pop_list_views(main_window())
        self.close()
        
    def okay(self):
        file_name=self.lineEdit_display_name.text()
        motor={
        'use_thermal': self.checkBox_use_thermal.isChecked(),
        'name': self.lineEdit_display_name.text(),
        'mass': self.doubleSpinBox_mass.value(),
        'Kv': self.doubleSpinBox_Kv.value(),
        'I0': self.doubleSpinBox_I0.value(),
        'Imax': self.doubleSpinBox_Imax.value(),
        'R': self.doubleSpinBox_R.value(),
        'Tref': self.doubleSpinBox_Tref.value(),
        'k1': self.doubleSpinBox_k1.value(),
        'k2': self.doubleSpinBox_k2.value(),
        'k3': self.doubleSpinBox_k3.value(),
        'Rthwm': self.doubleSpinBox_Rthwm.value(),
        'Rthma': self.doubleSpinBox_Rthma.value(),
        'Rthws': self.doubleSpinBox_Rthws.value(),
        'Rthsa': self.doubleSpinBox_Rthsa.value(),
        'alphar': self.doubleSpinBox_alphar.value(),
        'alpham': self.doubleSpinBox_alpham.value(),
        }
        
        
        try:   
            original_file_name=self.motor_edited['name']
            if motor==self.motor_edited:
                self.close()
                return
                
            QtWidgets.QMessageBox.warning(self, 'Motor edited',
                                            "Motor has been changed, its power surfaces must now be recalculated.",
                                            QtWidgets.QMessageBox.Ok)
            for root, dirs, files in os.walk('./power_surface_precalcs/'):
                for file in files:
                    if file.endswith(".pk") and motor['name'] in file:
                        os.remove('./power_surface_precalcs/' + file)
            try:
                os.remove('./motors/' + original_file_name + '.dat')
            except OSError:
                pass
        except AttributeError:
            pass
        
        data_dealings.create_dict_file('./motors/' + file_name + '.dat', motor, "// File defining motor constants and values to be read into dictionary")
        main_window.pop_list_views(main_window())
        
        self.close()
    
    def cancel(self):
        self.close()
        
class propeller_dialog_window(QtWidgets.QDialog):
    def __init__(self, propeller,prop_funcs, atmosphere, battery):
        super(propeller_dialog_window, self).__init__()
        uic.loadUi("./UI/dialog_prop_functions.ui", self)
        self.propeller=propeller
        self.prop_funcs=prop_funcs
        self.atmosphere=atmosphere
        self.battery=battery
        self.pushButton_calculate.clicked.connect(lambda: self.calculate(propeller, prop_funcs))
#        self.buttonBox.accepted.connect(lambda: edit_plane_window.okay(self))
        
    def populate(self, propeller):
        self.label_prop_name.setText(propeller['name'])
        self.label_diameter.setText("Diameter: {:0.2f} mm".format(propeller['diameter']*1000))
        self.doubleSpinBox_V.setValue(self.battery['V'])
        
    def calculate(self, propeller, prop_funcs):
        
        U=self.doubleSpinBox_U.value()
        Th=self.doubleSpinBox_Th.value()
        V=self.doubleSpinBox_V.value()
        
        if Th==0:
            QtWidgets.QMessageBox.warning(self, 'Thrust is 0',
                                            "Please input a non-zero propeller thrust",
                                            QtWidgets.QMessageBox.Ok)
            return
        
        atmosphere=self.atmosphere
        rho=atmosphere['rho']
        D=propeller['diameter']
#        Jrange=np.linspace(0,prop_funcs['bounds'][3],100)
#        C_T=prop_funcs['C_T_interp'](Jrange,Rerange)
#        U=prop_funcs['bounds'][3]*omega_search_max/np.pi*0.5*D
        
        omega_search=np.sqrt(Th/(0.05*rho*D**4))
        n=0
        nmax=10
        #This is an absolutely TERRIBLE way of doing this but I didn't have time to do it neatly
        while True:
            n+=1
            if n>nmax: #Stop if it doesn't converge
                QtWidgets.QMessageBox.warning(self, 'Unable to calculate operating point',
                                    "Check that your inputs are reasonable.",
                                    QtWidgets.QMessageBox.Ok)
                break
                return

            omega_search_min=0.5*omega_search-10
            omega_search_min=max([0,omega_search_min])
            omega_search_max=2*omega_search
            omega_range=np.linspace(omega_search_min,omega_search_max,100)
            
            Th_range=consumption_functions.propT(omega_range,U, prop_funcs['C_T_interp'], prop_funcs['bounds'], propeller,self.atmosphere)
#            print(Th_range)
            omega_int=interp1d(Th_range, omega_range, kind='linear')
            try:
                omega=float(omega_int(Th))
                break
            except ValueError as e:
#                print(e)
                if "A value in x_new is above the interpolation range." in str(e):
                    omega_search=omega_search_max
                elif "A value in x_new is below the interpolation range." in str(e):
                    omega_search=omega_search_min
#                break
        Q=consumption_functions.propQ(omega,U, prop_funcs['C_Q_interp'], prop_funcs['bounds'], propeller,self.atmosphere)[0]
#        print(Q)
        self.label_RPM.setText("Propeller RPM: {:0.2f}".format(omega*30/np.pi))
        self.label_Q.setText("Propeller torque: {:0.2f} mN m".format(Q*1000))
        
        Kv_min=omega*30/np.pi/V
        KT=V/omega
        Imax=Q/KT
        
        self.label_Kv.setText("Minimal Motor Kv: {:0.0f}".format(Kv_min))
        self.label_Imax.setText("Motor minimal current rating: {:0.2f}".format(Imax))
        
class kv_dialog_window(QtWidgets.QDialog):
    def __init__(self, atmosphere, battery):
        super(kv_dialog_window, self).__init__()
        uic.loadUi("./UI/dialog_kv_chooser.ui", self)
        self.atmosphere=atmosphere
        self.battery=battery
#        self.pushButton_calculate.clicked.connect(lambda: self.calculate(propeller, prop_funcs))
#        self.buttonBox.accepted.connect(lambda: edit_plane_window.okay(self))
        self.doubleSpinBox_mass.valueChanged.connect(self.refresh_thrust)
        self.spinBox_no_motors.valueChanged.connect(self.refresh_thrust)
        self.doubleSpinBox_Th.valueChanged.connect(self.refresh_mass)
        
        self.radioButton_manouv_low.toggled.connect(self.refresh_thrust)
        self.radioButton_manouv_med.toggled.connect(self.refresh_thrust)
        self.radioButton_manouv_high.toggled.connect(self.refresh_thrust)
        self.radioButton_manouv_custom.toggled.connect(self.refresh_thrust)
        
        self.low_factor=1.25
        self.med_factor=1.6
        self.high_factor=2.1

    def refresh_thrust(self):
        if self.doubleSpinBox_mass.value()==0:
            QtWidgets.QMessageBox.warning(self, 'Mass is 0',
                                            "Please input a non-zero craft mass",
                                            QtWidgets.QMessageBox.Ok)
            return
        Th_g=self.doubleSpinBox_mass.value()/self.spinBox_no_motors.value()*1000
        self.doubleSpinBox_Th.setValue(Th_g)
        
        if self.radioButton_manouv_low.isChecked():
            self.doubleSpinBox_Th_max.setValue(Th_g*self.low_factor)
            self.doubleSpinBox_Th_max.setEnabled(False)
        elif self.radioButton_manouv_med.isChecked():
            self.doubleSpinBox_Th_max.setValue(Th_g*self.med_factor)
            self.doubleSpinBox_Th_max.setEnabled(False)
        elif self.radioButton_manouv_high.isChecked():
            self.doubleSpinBox_Th_max.setValue(Th_g*self.high_factor)
            self.doubleSpinBox_Th_max.setEnabled(False)
        elif self.radioButton_manouv_custom.isChecked():
            self.doubleSpinBox_Th_max.setEnabled(True)
            

    def refresh_mass(self):
        if self.doubleSpinBox_Th.value()==0:
            QtWidgets.QMessageBox.warning(self, 'Thrust is 0',
                                            "Please input a non-zero thrust",
                                            QtWidgets.QMessageBox.Ok)
            return
        mass_kg=self.doubleSpinBox_Th.value()/1000*self.spinBox_no_motors.value()
        self.doubleSpinBox_mass.setValue(mass_kg)
        Th_g=self.doubleSpinBox_Th.value()
        
        if self.radioButton_manouv_low.isChecked():
            self.doubleSpinBox_Th_max.setValue(Th_g*self.low_factor)
        elif self.radioButton_manouv_med.isChecked():
            self.doubleSpinBox_Th_max.setValue(Th_g*self.med_factor)
        elif self.radioButton_manouv_high.isChecked():
            self.doubleSpinBox_Th_max.setValue(Th_g*self.high_factor)

    def populate(self):

        self.doubleSpinBox_V.setValue(self.battery['V'])
        
    def calculate(self):
        
        Th=self.doubleSpinBox_Th.value()*9.8*1000
        
        if Th==0:
            QtWidgets.QMessageBox.warning(self, 'Thrust is 0',
                                            "Please input a non-zero propeller thrust",
                                            QtWidgets.QMessageBox.Ok)
            return
        
        atmosphere=self.atmosphere
        rho=atmosphere['rho']
        
        
if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = main_window()
    window.show()
    sys.exit(app.exec_())