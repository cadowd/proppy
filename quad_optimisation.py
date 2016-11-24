# -*- coding: utf-8 -*-
"""
Created on Nov 22 14:51:57 2016

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

from PyQt4 import QtCore, QtGui, uic
import data_dealings
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
import numpy as np
from matplotlib.path import Path
import matplotlib.patches as patches
import consumption_functions
import scenarios
import main
from scipy.interpolate import interp1d
from operator import itemgetter
from matplotlib import cm
import itertools
import time

def do_calc(self, motor_dat, prop_dat, thrust):
    """
    Calculates static performance of combination.
    """
    C=np.linspace(0,1.0,100)
    motor=data_dealings.read_dict_file('./motors/' + motor_dat + '.dat')
    propeller=main.get_prop_dict(self, prop_dat)
    
    Imax, rpm_static, Tmax=scenarios.static_max_current(C,self.atmosphere,self.battery,motor,propeller)
    
#    print(Tmax)
    I_interp=interp1d(Tmax, Imax, kind='linear')
#    print(thrust)
#    print(np.nanmax(Tmax))
    try:
        I_thrust=I_interp(thrust)
    except ValueError:
#        print(e)
        I_thrust=float('nan')
    Pel_nom=self.battery['V']*I_thrust
    Pel=self.battery['V']*Imax
    result={
    'I': Imax,
    'rpm': rpm_static,
    'T': Tmax,
    'Imax': np.nanmax(Imax),
    'Tmax': np.nanmax(Tmax),
    'Pel_nom': Pel_nom,
    'Pel': Pel
    }
    
    
    return result


class quad_optimise_window(QtGui.QDialog):
    def __init__(self, propellers, motors, options):
        super(quad_optimise_window, self).__init__()
        uic.loadUi("./UI/quad_optimisation_window.ui", self)
        self.atmosphere=options['atmosphere']
        self.battery=options['battery']
        self.doubleSpinBox_voltage.setValue(self.battery['V'])
        self.Urange=np.linspace(0.5,30,100)
        fig1 = Figure()
        ax1f1 = fig1.add_subplot(111)
        self.canvas = FigureCanvas(fig1)
        self.layout_plot.addWidget(self.canvas)
        self.toolbar = NavigationToolbar(self.canvas, 
                self, coordinates=True)
#        self.addToolBar(self.toolbar)
        self.layout_plot.addWidget(self.toolbar)
        
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
        
        self.pushButton_optimise.clicked.connect(lambda: self.get_optimal(ax1f1))
        
    def get_optimal(self, axes):
        """
        Function to analyse and compare ALL the selected combinations.
        """
        axes.clear()
        motors_dat=main.return_checked_values(self.listView_motors)
        propellers_dat=main.return_checked_values(self.listView_props)
        
        self.battery['V']=self.doubleSpinBox_voltage.value()
        no_motors=self.spinBox_no_motors.value()
        number=self.spinBox_number.value()
        mass=self.doubleSpinBox_mass.value()
        thrust_optimise=mass*9.8/no_motors
        
        battery=self.battery
        atmosphere=self.atmosphere
        
        min_acc=self.doubleSpinBox_min_acc.value()
        max_I_perc=self.doubleSpinBox_max_I.value()/100
        
        if mass==0:
            QtGui.QMessageBox.warning(self, 'No craft mass',
                                            "Input the mass of the craft.",
                                            QtGui.QMessageBox.Ok)
            return
        
        no_to_plot=0
        calc_list=[]
        for motor_dat in motors_dat:
            for prop_dat in propellers_dat:
#                line_label=motor_dat + ' ' + prop_dat
                no_to_plot +=1
                calc_dict={
                'motor':motor_dat,
                'prop':prop_dat,
                'label': motor_dat + ' ' + prop_dat}
                calc_list.append(calc_dict)
                
        if len(calc_list)==0:
            return
        self.prog_widget = QtGui.QProgressDialog("Optimising selection", "Stop the madness!",0,no_to_plot)
        self.prog_widget.setWindowModality(QtCore.Qt.WindowModal)
#        self.prog_widget.setMinimumDuration(0)
        self.prog_widget.show()

        Capacity=self.battery['capacity']*3.6*self.battery['V']
        
        # Setting the value on every run to 0
        self.prog_widget.setValue(0)
        result_list=[]
        for calc_dict in calc_list:
            QtGui.QApplication.instance().processEvents()
            if self.prog_widget.wasCanceled():
                QtGui.QMessageBox.information(self, 'Optimisation cancelled',
                                        "The optimisation has been cancelled.",
                                        QtGui.QMessageBox.Ok)
                return

            result = do_calc(self, calc_dict['motor'], calc_dict['prop'], thrust_optimise)
            self.prog_widget.setValue(self.prog_widget.value()+1)
            motor_dict=data_dealings.read_dict_file('./motors/' + calc_dict['motor'] + '.dat')
#            propeller=main.get_prop_dict(self, prop_dat)
            print(result['Pel_nom'])
            if np.isnan(result['Pel_nom']):
                continue
                print('No power')
            if result['Imax']>max_I_perc*motor_dict['Imax']:
                print('overcurrent')
                continue
            max_acc=result['Tmax']*no_motors/mass - mass*9.8
            if max_acc<min_acc:
                print('Low acc ' + str(max_acc) + "m/s/s")
                continue
            result['label']=calc_dict['label']
            result['max_time']=int(Capacity/result['Pel_nom'])
            result['max_acc']=max_acc
            result_list.append(result)
        
        if len(result_list)<1:
            QtGui.QMessageBox.warning(self, 'No suitable combinations found',
                                            "No combinations were found that match your constraints. Consider loosening the constraints or selecting different components.",
                                            QtGui.QMessageBox.Ok)
            return
#            self.single_plot(result, axes)
        axes.set_xlabel(r"System thrust per motor[N]")
        axes.set_title("Power consumption")
        axes.set_ylabel(r"Power consumed [W]")

        top_combos=sorted(result_list, key=itemgetter('max_time'), reverse=True)[:number]
        
        combo_list = QtGui.QStandardItemModel(self.listView_results)
        for combo in top_combos:
            label="{} - Max hover time: {} - Max upwards acceleration: {:0.2f}m/s/s ".format(combo['label'], 
            time.strftime('%Hh%Mm%Ss', time.gmtime(combo['max_time'])), combo['max_acc'])
            print(label)
            item = QtGui.QStandardItem(label)
            item.setCheckable(True)
            item.setCheckState(QtCore.Qt.Checked)
            item.combo=combo
            combo_list.appendRow(item)
        self.listView_results.setModel(combo_list)
        
        
        idx = np.linspace(0, 1, number)
        cmap = cm.get_cmap('Accent')
        colors = itertools.cycle(cmap(idx))
        axes.axvline(thrust_optimise)
        for result in top_combos:
            c=next(colors)
            axes.plot(result['T'],result['Pel'], label=result['label'], color=c)
        axes.legend(loc=4)
        self.canvas.draw()