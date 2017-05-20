# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 15:48:20 2017

@author: lis-15-15
"""

import sys
from PyQt5 import QtCore, QtGui, uic, QtWidgets, Qt
#import UI.qtRangeSlider as qtRangeSlider
import UI.RangeSlider as RangeSlider

import os
from matplotlib.figure import Figure
import itertools
import numpy as np
import scipy as sp
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
    
import pickle
import numbers
from scipy.interpolate import interp1d

import data_dealings
 
qtCreatorFile = "./UI/motor_selection_window.ui" # Enter file here.
 
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

def get_csv_names(folder):
    #Get names and locations of all suitable data files
    file_names=[]
    for file in os.listdir(folder):
        if file.endswith(".csv"):
#            data_file=os.path.join(folder, file)
            file_name=file.split("/")[-1]
            file_name=file_name.rsplit(".",1)[0]
            file_names.append(file_name)
                
    return file_names

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

def return_checked_values(listView):
    checked_items=[]
    
    item_model = listView.model()
    for row in range(item_model.rowCount()):
        item = item_model.item(row)
        if item.checkState() == QtCore.Qt.Checked:
            checked_items.append(item.text())
    return checked_items

class main_window(Ui_MainWindow, QtBaseClass):
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)

        self.setupUi(self)
        
#        fig1 = Figure()
#        ax1f1 = fig1.add_subplot(111)
#        self.canvas = FigureCanvas(fig1)
#        self.layout_plot.addWidget(self.canvas)
#        self.toolbar = NavigationToolbar(self.canvas, 
#                self, coordinates=True)
#        self.addToolBar(self.toolbar)
#        self.axes=ax1f1
#        self.layout_plot.addWidget(self.toolbar)
        self.plot_colors = itertools.cycle(['k', "r", "b", "g", "y", "c","m"])
        self.current_plots = []
        self.current_data = []
        self.current_data_files = []
        self.file_folder="."
        
        self.get_full_motor_list()
        self.header_list=["name", "Kv", "mass", "max amps", "supplier"]
        self.motor_list_model = QtGui.QStandardItemModel(self.tableView_motor_list)
        self.motor_list_model.setHorizontalHeaderLabels(self.header_list)
        
#        for header in self.header_list:
#            motor_list.setHeaderData(header, Qt.Horizontal, QVariant(header))
        
        self.maxA_global=0
        self.minA_global=0
        self.maxKv_global=0
        self.minKv_global=0
        self.maxMass_global=0
        self.minMass_global=0
        self.pop_table()
        
        self.rangeSlider_Kv = RangeSlider.RangeSlider(Qt.Qt.Horizontal)
        layout = QtWidgets.QGridLayout(self.rangeSliderWidget_Kv)
#        layout.setMargin(1)
        layout.addWidget(self.rangeSlider_Kv)
        self.rangeSlider_Kv.setMinimum(self.minKv_global)
        self.rangeSlider_Kv.setMaximum(self.maxKv_global)
        self.rangeSlider_Kv.setLow(self.minKv_global)
        self.rangeSlider_Kv.setHigh(self.maxKv_global)
        self.rangeSlider_Kv.setTickPosition(QtWidgets.QSlider.TicksBelow)
        
        self.rangeSlider_mass = RangeSlider.RangeSlider(Qt.Qt.Horizontal)
        layout = QtWidgets.QGridLayout(self.rangeSliderWidget_mass)
#        layout.setMargin(1)
        layout.addWidget(self.rangeSlider_mass)
        self.rangeSlider_mass.setMinimum(self.minMass_global)
        self.rangeSlider_mass.setMaximum(self.maxMass_global)
        self.rangeSlider_mass.setLow(self.minMass_global)
        self.rangeSlider_mass.setHigh(self.maxMass_global)
        self.rangeSlider_mass.setTickPosition(QtWidgets.QSlider.TicksBelow)
        
        self.rangeSlider_maxA = RangeSlider.RangeSlider(Qt.Qt.Horizontal)
        layout = QtWidgets.QGridLayout(self.rangeSliderWidget_maxA)
#        layout.setMargin(1)
        layout.addWidget(self.rangeSlider_maxA)
        self.rangeSlider_maxA.setMinimum(0)
        self.rangeSlider_maxA.setMaximum(self.maxA_global)
        self.rangeSlider_maxA.setLow(0)
        self.rangeSlider_maxA.setHigh(self.maxA_global)
        self.rangeSlider_maxA.setTickPosition(QtWidgets.QSlider.TicksBelow)

        self.label_kv_max.setText("{:.0f}".format(self.maxKv_global))
        self.label_kv_min.setText("{:.0f}".format(self.minKv_global))
        self.label_mass_min.setText("{:.0f}".format(self.minMass_global))
        self.label_mass_max.setText("{:.0f}".format(self.maxMass_global))
        self.label_amps_min.setText("{:.0f}".format(self.minA_global))
        self.label_amps_max.setText("{:.0f}".format(self.maxA_global))
        
        
#        #Button functions and so on
        self.rangeSlider_Kv.sliderMoved.connect(self.kv_sliderMoved)
        self.rangeSlider_mass.sliderMoved.connect(self.mass_sliderMoved)
        self.rangeSlider_maxA.sliderMoved.connect(self.maxA_sliderMoved)
#        self.pushButton_plot.clicked.connect(lambda: self.plot(ax1f1))
#        self.pushButton_clear_plot.clicked.connect(lambda: self.clear_plot(ax1f1))
#        
#        self.radioButton_V_Pel.toggled.connect(lambda: self.change_plot(ax1f1))
#        self.radioButton_V_eta.toggled.connect(lambda: self.change_plot(ax1f1))
#        self.radioButton_V_eta_prop.toggled.connect(lambda: self.change_plot(ax1f1))
#        self.radioButton_V_eta_mot.toggled.connect(lambda: self.change_plot(ax1f1))
#        
#        
#        self.pushButton_motor_edit.clicked.connect(self.edit_motor)
#        self.pushButton_motor_create.clicked.connect(self.new_motor)
#        
#        self.pushButton_plane_edit.clicked.connect(self.edit_plane)
#        self.pushButton_plane_create.clicked.connect(self.create_plane)
#        
#        self.pushButton_prop_filter.clicked.connect(self.filter_props)
#        
#        
#        self.pushButton_optimise.clicked.connect(self.optimise_window)
#        self.pushButton_quad_optimise.clicked.connect(self.quad_optimise_window)
#        self.pushButton_get_risks.clicked.connect(self.get_risks)
#        self.pushButton_prop_dets.clicked.connect(self.propeller_functions)
#        
#        self.opt_settings=pickle.load( open( 'opt_settings.pk', "rb" ) )

        
        self.battery={
        'V':14.8, #battery voltage
        'capacity':4800
        }
        
        self.atmosphere={
        'rho':1.20, #
        'mu':1.8e-5, #dynamic viscosity
        'Ta': 25, #Ambient temperature
        }
        

        
    def kv_sliderMoved(self):
        lower=self.rangeSlider_Kv.low()
        upper=self.rangeSlider_Kv.high()
        label=str(lower)+" - "+str(upper)
        self.label_kv_range.setText(label)
        
    def mass_sliderMoved(self):
        lower=self.rangeSlider_mass.low()
        upper=self.rangeSlider_mass.high()
        label=str(lower)+" - "+str(upper)
        self.label_mass_range.setText(label)
        
    def maxA_sliderMoved(self):
        lower=self.rangeSlider_maxA.low()
        upper=self.rangeSlider_maxA.high()
        label=str(lower)+" - "+str(upper)
        self.label_amps_range.setText(label)
        
    def get_full_motor_list(self):
        motor_list=[]
        header_list=[]
        for name in get_dat_names('./motors/'):
            motor_dat=data_dealings.read_dict_file('./motors/' + name + '.dat')
            motor_list.append(motor_dat)
            if len(motor_dat.keys())>len(header_list):
                header_list=[]
                for key in motor_dat.keys():
                    header_list.append(key)
        self.motor_list=motor_list
        self.header_list=header_list
        
    def pop_table(self):
        """
        Function to populate lists with the files found in the folders.
        """
        
        for motor in self.motor_list:
#            motor_dat=data_dealings.read_dict_file('./motors/' + name + '.dat')
            #Find global maxs and mins
            if motor["Imax"]>self.maxA_global:
                self.maxA_global=motor["Imax"]
            if motor["Imax"]<self.minA_global:
                self.minA_global=motor["Imax"]
            if motor["mass"]>self.maxMass_global:
                self.maxMass_global=motor["mass"]
            if motor["mass"]<self.minMass_global:
                self.minMass_global=motor["mass"]
            if motor["Kv"]>self.maxKv_global:
                self.maxKv_global=motor["Kv"]
            if motor["Kv"]<self.minKv_global:
                self.minKv_global=motor["Kv"]
                
            item_list=[]
            for header in self.header_list:
                try:
                    print(motor[header])
                    item = QtGui.QStandardItem(str(motor[header]))
                except KeyError:
                    if header == "max amps":
                        item = QtGui.QStandardItem(str(motor["Imax"]))
                    else:
                        item = QtGui.QStandardItem("")
                if header == "name":
                    item.setCheckable(True)
                    item.data=motor
#                item.s
                item_list.append(item)
#            item.data = motor_dat
#            print(item_list)
            self.motor_list_model.appendRow(item_list)
        
        self.motor_list_proxy=QtCore.QSortFilterProxyModel(self)
        self.motor_list_proxy.setSourceModel(self.motor_list_model)
        
        self.tableView_motor_list.setModel(self.motor_list_proxy)
#        self.tableView_motor_list.setModel(self.motor_list_model)

    def filter_table(self):
        self.motor_list_proxy.lessThan()
                

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = main_window()
    window.show()
    sys.exit(app.exec_())