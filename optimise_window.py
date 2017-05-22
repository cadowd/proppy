# -*- coding: utf-8 -*-
"""
Created on Nov 23 11:13:02 2016

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
import time
from PyQt5 import QtCore, QtGui, uic, QtWidgets
import os
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import itertools
import numpy as np
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
    
import consumption_functions as optimised_consumption
import data_dealings
import main
import comparison_functions
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib import path, cm, rcParams, collections

qtCreatorFile = "./UI/main_window.ui" # Enter file here.
 
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

#dialog_edit_motor = uic.loadUiType("./UI/dialog_edit_motor.ui")

    
class optimise_window(QtWidgets.QDialog):
    def __init__(self, plane, propellers, motors, options):
        super(optimise_window, self).__init__()
        uic.loadUi("./UI/optimisation_window.ui", self)
        self.plane_label.setText(plane)
        plane_font=QtGui.QFont("Arial", 14, QtGui.QFont.Bold)
        self.plane_label.setFont(plane_font)
#
#        self.setupUi(self)
        self.battery=options['battery']
        self.atmosphere=options['atmosphere']
        self.optimisation=options['optimisation']
        plane_dict=data_dealings.read_dict_file('./planes/' + plane + '.dat')
        self.plane=plane_dict
        
        fig1 = Figure()
        ax1f1 = fig1.add_subplot(111)
        self.canvas = FigureCanvas(fig1)
        self.layout_plot.addWidget(self.canvas)
        self.toolbar = NavigationToolbar(self.canvas, 
                self, coordinates=True)
#        self.addToolBar(self.toolbar)
        self.layout_plot.addWidget(self.toolbar)

#        
        prop_list = QtGui.QStandardItemModel(self.listView_props)
        for name in propellers:
            item = QtGui.QStandardItem(name)
            item.setCheckable(True)
            item.setCheckState(QtCore.Qt.Checked)
            prop_list.appendRow(item)
        self.listView_props.setModel(prop_list)
        
        motor_list = QtGui.QStandardItemModel(self.listView_props)
        for name in motors:
            item = QtGui.QStandardItem(name)
            item.setCheckable(True)
            item.setCheckState(QtCore.Qt.Checked)
            motor_list.appendRow(item)
        self.listView_motors.setModel(motor_list)
       
        self.radioButton_optimal_speed.toggled.connect(self.refresh_window)
        self.radioButton_fixed_speed.toggled.connect(self.refresh_window)
        
        self.checkBox_static_Th_const.clicked.connect(self.refresh_window)
        self.checkBox_max_current_const.clicked.connect(self.refresh_window)
        self.checkBox_top_speed_const.clicked.connect(self.refresh_window)
        self.checkBox_motor_max_mass_const.clicked.connect(self.refresh_window)
        
        self.pushButton_cancel.clicked.connect(self.cancel)
        self.pushButton_get_optimal.clicked.connect(lambda: self.get_optimal(plane, ax1f1))
        
        self.pushButton_plot_performance.clicked.connect(lambda: self.plot_perf(ax1f1))
        self.pushButton_plot_efficiency.clicked.connect(lambda: self.plot_eff(ax1f1))
        
        


    def cancel(self):
        self.close()
        
    def refresh_window(self):
        """
        Function refresh the window when a checkbox is changed.
        """
        if self.radioButton_optimal_speed.isChecked():
            self.doubleSpinBox_flight_speed.setEnabled(False)
        else:
            self.doubleSpinBox_flight_speed.setEnabled(True)
        
        if self.checkBox_static_Th_const.isChecked():
            self.doubleSpinBox_static_Th.setEnabled(True)
        else:
            self.doubleSpinBox_static_Th.setEnabled(False)
        
        if self.checkBox_max_current_const.isChecked():
            self.doubleSpinBox_max_current.setEnabled(True)
        else:
            self.doubleSpinBox_max_current.setEnabled(False)
        
        if self.checkBox_top_speed_const.isChecked():
            self.doubleSpinBox_top_speed.setEnabled(True)
        else:
            self.doubleSpinBox_top_speed.setEnabled(False)
        
        if self.checkBox_motor_max_mass_const.isChecked():
            self.doubleSpinBox_motor_max_mass.setEnabled(True)
        else:
            self.doubleSpinBox_motor_max_mass.setEnabled(False)
            
    def get_optimal(self, plane_dat, axes):
        """
        Function to analyse and compare ALL the selected combinations.
        """
        motors_dat=main.return_checked_values(self.listView_motors)
        propellers_dat=main.return_checked_values(self.listView_props)
        
        Urange= np.linspace(self.optimisation['Umin'], self.optimisation['Umax'], self.optimisation['resolution'])
        Trange = np.linspace(self.optimisation['Tmin'], self.optimisation['Tmax'], self.optimisation['resolution'])
        folder='./power_surface_precalcs/'
        motors=[]
        propellers=[]
        
        battery=self.battery
        atmosphere=self.atmosphere
        
        for motor_dat in motors_dat:
            motors.append(data_dealings.read_dict_file('./motors/' + motor_dat + '.dat'))
        
        for propeller_dat in propellers_dat:
            propellers.append(main.get_prop_dict(self, propeller_dat))

        plane=self.plane
        
        
        
        combo_list=comparison_functions.check_Psurfs(motors,propellers,Urange, Trange, folder, battery, atmosphere)
        if len(combo_list)>0:
            progress_widgit=QtWidgets.QProgressDialog("Generating missing power surfaces", "Stop the madness!",0,len(combo_list))
            progress_widgit.setWindowModality(QtCore.Qt.WindowModal)
            
            progress_widgit.show()
            progress_widgit.setValue(0)
            for combo_dict in combo_list:
                comparison_functions.gen_single_Psurf(combo_dict,Urange, Trange, folder, battery, atmosphere)
                progress_widgit.setValue(progress_widgit.value()+1)
                if progress_widgit.wasCanceled():
                    QtWidgets.QMessageBox.information(self, 'Optimisation cancelled',
                                            "The optimisation has been cancelled.",
                                            QtWidgets.QMessageBox.Ok)
                    return
        
        
        conditions={
        'maximise_distance': True,
        'fixed_speed': False, #10.0, #False or number
        'top_speed': False,
        'static_thrust': False,
        'allowable_overload': False,
        'motor_mass': False,
        }
        
        #Don't worry about highlighting the best regions if there's only 1 selection.
        if len(motors)*len(propellers)>1:
            best_regions=comparison_functions.calc_best_regions(motors, propellers, folder, conditions)
            self.plot_regions( axes, best_regions, plane)
        
        
        if self.radioButton_maximise_time.isChecked():
            conditions['maximise_distance']=False
        if self.radioButton_fixed_speed.isChecked():
            conditions['fixed_speed']=self.doubleSpinBox_flight_speed.value()
        if self.checkBox_static_Th_const.isChecked():
            conditions['static_thrust']=self.doubleSpinBox_static_Th.value()
        if self.checkBox_max_current_const.isChecked():
            conditions['allowable_overload']=self.doubleSpinBox_max_current.value()/100
        if self.checkBox_top_speed_const.isChecked():
            conditions['top_speed']=self.doubleSpinBox_top_speed.value()
        if self.checkBox_motor_max_mass_const.isChecked():
            conditions['motor_mass']=self.doubleSpinBox_motor_max_mass.value()
        
        no_tops=self.spinBox_no_best.value()
        top_combos, report=comparison_functions.top_choices(plane, motors, propellers, atmosphere, folder, conditions, no_tops)
        print(len(top_combos))
        
        QtWidgets.QMessageBox.information(self, 'Analysis complete.',
                                """Combinations analysed: {} \n
Take off thrust too low: {} \n
Maximum current too high: {} \n
Over motor mass limit: {} \n
Unable to produce enough thrust for flight: {} \n
Suitable: {} \n""".format(report['analysed'], report['low_thrust'], report['over_current'], report['over_mass'],report['too_weak'],report['suitable']),
                                QtWidgets.QMessageBox.Ok)
        
        combo_list = QtGui.QStandardItemModel(self.listView_results)
        for combo in top_combos:
            label="{} - Max flight time: {} - Max distance: {:0.2f}m - Cruise speed: {:0.2f}m/s".format(combo['combo'], 
            time.strftime('%Hh%Mm%Ss', time.gmtime(combo['max_time'])), combo['max_distance'], combo['max_time_speed'])
            print(label)
            item = QtGui.QStandardItem(label)
            item.setCheckable(True)
            item.setCheckState(QtCore.Qt.Checked)
            item.combo=combo
            combo_list.appendRow(item)
        self.listView_results.setModel(combo_list)
        
        self.optmisation_results=top_combos

        
        
        
    def plot_regions(self, axes, best_regions, plane):
        """
        Function to plot the best efficiency regions, once they've been calculated.
        """
        axes.autoscale(True)
#        colors = itertools.cycle(["r", "b", "g", "y", "c","m"])
        idx = np.linspace(0, 1, len(best_regions))
        cmap = cm.get_cmap('Accent')
        colors = itertools.cycle(cmap(idx))
        
        legend_handles=[]
        axes.clear()
        axes.axis([self.optimisation['Umin'], self.optimisation['Umax'], self.optimisation['Tmin'], self.optimisation['Tmax']])

        
        battery=self.battery
        atmosphere=self.atmosphere
        
#        fig, ax = plt.subplots()
        
#        axes.set_color_cycle(jet(idx))
        
        for patch in best_regions:
            zone=patch['region']
            combo_colour=next(colors)
            area=0
            for i,p in enumerate(zone):
                axes.add_collection(collections.PathCollection(zone, color=combo_colour))
                bbx=p.get_extents()
                length=bbx.x1 - bbx.x0
                height=bbx.y1 - bbx.y0
                area=area+length*height

            if area>0.1:
                combo_patch = mpatches.Patch(color=combo_colour, label=patch['combo_name'])
                legend_handles.append(combo_patch)
                    
        
#        plt.scatter(points[:,0], points[:,1])
        
        axes.set_title(r"Most efficient combination")
        axes.set_ylabel(r"Thrust produced [N]")
        axes.set_xlabel(r"Air speed (in steady flight) [m/s]")
        
        axes.autoscale(False)
        Uplot=np.linspace(self.optimisation['Umin']+0.1, self.optimisation['Umax'], 100)
        dragP=optimised_consumption.dragFunc(Uplot, plane, atmosphere)
        self.plane_curve=[Uplot, dragP]
        axes.plot(Uplot,dragP, color='k', label='Plane')
        plane_line = mlines.Line2D([], [], color='k', label='Plane drag')
        
        legend_handles.append(plane_line)
        axes.legend(handles=legend_handles)
        self.canvas.draw()
        
    def plot_perf(self, axes):
        """
        Function to plot performance of currently selected result.
        """
        selected_combos=[]
        axes.clear()
    
        item_model = self.listView_results.model()
        for row in self.listView_results.selectedIndexes():
            item = item_model.itemFromIndex(row)
            selected_combos.append(item.combo)
            
        combo=selected_combos[0]     
        plane=self.plane
        atmosphere=self.atmosphere
        
#        plt.suptitle(image_name, fontsize=14, fontweight='bold')
        axes.set_title(r"Power consumption of system")
        axes.set_ylabel(r"Thrust produced [N]")
        axes.set_xlabel(r"Air speed (in steady flight) [m/s]")
        Pmax=combo['motor']['Imax']*self.battery['V']
        levels = np.linspace(20,Pmax,10) #[20,40,60,80,100,120,140,160,180,200,220,240,260,280]
        colormap=cm.get_cmap('coolwarm')
        CS_filled=axes.contourf(combo['Umesh'], combo['Tmesh'], combo['Pmesh'],levels, cmap=colormap, extend='both')
        if hasattr(self, 'cbar'):
#            image.set_data(CS_filled) #use this if you use new array
            self.cbar = plt.colorbar(CS_filled, cax=self.cbar.ax)
            self.cbar.ax.set_ylabel('Power [W]')
        else:
            self.cbar = plt.colorbar(CS_filled, ax=axes)
            self.cbar.ax.set_ylabel('Power [W]')

        CS=axes.contour(combo['Umesh'], combo['Tmesh'], combo['Pmesh'],levels, colors=('k',), linewidths=(1,))
        axes.clabel(CS, fmt='%2.1f', colors='k', fontsize=14)
        
        axes.contourf(combo['Umesh'], combo['Tmesh'], combo['Pmesh'],[Pmax,Pmax*1.2], hatches=['//'], colors='none', extend='max',alpha=0.0)
        #Find the max because autoscale can't handle the abundunce of nana in these arrays.
        id_valid=(np.array(combo['Pmesh']) == np.array(combo['Pmesh']))
        T=combo['Tmesh']
        T[~id_valid]=float('nan')
        axes.axis([0, np.nanmax(combo['Umesh']), 0, np.nanmax(T)])
#        axes.autoscale(True)
        if self.radioButton_fixed_speed.isChecked():
            axes.scatter(combo['max_time_speed'], combo['max_time_drag'], c='b', label='Cruise speed')
        else:
            axes.scatter(combo['max_distance_speed'], combo['max_distance_drag'], c='b', label='Max range conditions')
            axes.scatter(combo['max_time_speed'], combo['max_time_drag'], c='r', label='Max flight time conditions')
        try:
            axes.plot(self.plane_curve[0],self.plane_curve[1], color='k', label='Plane drag')
        except AttributeError:
            Uplot=np.linspace(self.optimisation['Umin']+0.1, self.optimisation['Umax'], 100)
            dragP=optimised_consumption.dragFunc(Uplot, plane, atmosphere)
            self.plane_curve=[Uplot, dragP]
            axes.plot(Uplot,dragP, color='k', label='Plane drag')
            plane_line = mlines.Line2D([], [], color='k', label='Plane drag')
        
        axes.legend(loc=1)
        self.canvas.draw()
        
    def plot_eff(self, axes):
#        try:
#            plt.delaxes(axes.axes[1])
#        except Exception as e:
#            print(e)
        selected_combos=[]
        axes.clear()
        
    
        item_model = self.listView_results.model()
        for row in self.listView_results.selectedIndexes():
            item = item_model.itemFromIndex(row)
            selected_combos.append(item.combo)
            
        combo=selected_combos[0]        
        plane=self.plane
        atmosphere=self.atmosphere

#            plt.suptitle(image_name, fontsize=14, fontweight='bold')
        axes.set_title(r"Total efficiency of combination")
        axes.set_ylabel(r"Thrust produced [N]")
        axes.set_xlabel(r"Air speed (in steady flight) [m/s]")
        levels = [0.3,0.35,0.4,0.45,0.5,0.55,0.6, 0.65, 0.70, 0.75]
        colormap=cm.get_cmap('RdYlGn')
        CS_filled=axes.contourf(combo['Umesh'], combo['Tmesh'], combo['eta_system'],levels, cmap=colormap, extend='both')

        if hasattr(self, 'cbar'):
            self.cbar = plt.colorbar(CS_filled, cax=self.cbar.ax)
            self.cbar.ax.set_ylabel('Total efficiency')
        else:
            self.cbar = plt.colorbar(CS_filled, ax=axes)
            self.cbar.ax.set_ylabel('Total efficiency')
        

        CS=plt.contour(combo['Umesh'], combo['Tmesh'], combo['eta_system'], levels, colors=('k',), linewidths=(1,))
        axes.clabel(CS, fmt='%2.1f', colors='k', fontsize=14)
        
        id_valid=(np.array(combo['Pmesh']) == np.array(combo['Pmesh']))
        T=combo['Tmesh']
        T[~id_valid]=float('nan')
        axes.axis([0, np.nanmax(combo['Umesh']), 0, np.nanmax(T)])
#        axes.plot(self.plane_curve[0],self.plane_curve[1], color='k', label='PLane')
        if self.radioButton_fixed_speed.isChecked():
            axes.scatter(combo['max_time_speed'], combo['max_time_drag'], c='b', label='Cruise speed')
        else:
            axes.scatter(combo['max_distance_speed'], combo['max_distance_drag'], c='b', label='Max range conditions')
            axes.scatter(combo['max_time_speed'], combo['max_time_drag'], c='r', label='Max flight time conditions')

        try:
            axes.plot(self.plane_curve[0],self.plane_curve[1], color='k', label='Plane drag')
        except AttributeError:
            Uplot=np.linspace(self.optimisation['Umin']+0.1, self.optimisation['Umax'], 100)
            dragP=optimised_consumption.dragFunc(Uplot, plane, atmosphere)
            self.plane_curve=[Uplot, dragP]
            axes.plot(Uplot,dragP, color='k', label='Plane drag')
            plane_line = mlines.Line2D([], [], color='k', label='Plane drag')

        axes.legend(loc=1)
        self.canvas.draw()
        
        