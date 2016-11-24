# -*- coding: utf-8 -*-
"""


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

class estimate_plane_window(QtGui.QDialog):
    def __init__(self, atmosphere):
        super(estimate_plane_window, self).__init__()
        uic.loadUi("./UI/estimator_dialogue.ui", self)
        self.atmosphere=atmosphere
        self.Urange=np.linspace(0.5,30,100)
        fig1 = Figure()
        shape = fig1.add_subplot(111)
        fig1.patch.set_visible(False)
        shape.set_autoscale_on(False)
#        shape.axis('off')
        self.canvas_shape = FigureCanvas(fig1)
        self.layout_shape.addWidget(self.canvas_shape)

        
        verts = [
        (-0.1, 0), # left, bottom
        (-0.1, 0.1), # left, top
        (0, 0.1), # right, top
        (0.1, 0.1), # right, bottom
        (0.1, 0), # ignored
        (0, 0), # ignored
        (0, 0), # ignored
        ]
    
        codes = [Path.MOVETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.CLOSEPOLY,
                 ]
        
        path = Path(verts, codes)
        
#        self.plane_patch = patches.PathPatch(path, facecolor='gray', lw=2)
        self.canvas_shape.draw()
        self.background = self.canvas_shape.copy_from_bbox(shape.bbox)
        
        self.plane_patch = patches.Polygon(verts, lw=2, facecolor='gray', closed ='True', animated=True)
        shape.axis([-1.5, 1.5, -1.5, 1.5])
#        shape.add_patch(self.plane_patch)
        shape.add_artist(self.plane_patch)
        self.xlim_old=1
#        self.canvas_shape.draw()
        
        
        fig2 = Figure()
        polar = fig2.add_subplot(111)
        self.canvas_polar = FigureCanvas(fig2)
        self.layout_polar.addWidget(self.canvas_polar)
        self.toolbar = NavigationToolbar(self.canvas_polar, 
                self, coordinates=True)
#        self.addToolBar(self.toolbar)
#        self.axes=ax1f1
        self.layout_polar.addWidget(self.toolbar)
        
        polar.set_xlabel(r"Airspeed in steady flight")
        polar.set_title("Drag vs air speed")
        polar.set_ylabel(r"Drag [N]")
        self.polar_line=polar.plot(self.Urange, self.Urange, label='Estimated plane')[0]
        polar.axis([0, 30, 0, 5])
        
        
        polar.set_xlabel(r"Air speed (in steady flight) [m/s]")
        
        self.rho=self.atmosphere['rho']
        self.payload_ratio_min=0.08 #You could probably do worse than this if you wanted, but its not that interesting
        self.payload_ratio_max=0.4 #Best I could find were airframes such as the Maja and the X8 skywalker
        
        structure_factor=self.verticalSlider_structure_factor.value()/50
        self.aero_factor=self.verticalSlider_aero_factor.value()/50
        self.theta_sweep=self.doubleSpinBox_theta_sweep.value()*np.pi/180
        self.b_ratio=self.horizontalSlider_span.value()/100
        self.payload_ratio=self.payload_ratio_min + structure_factor*(self.payload_ratio_max-self.payload_ratio_min)        
        self.Q=self.horizontalSlider_Q_factor.value()/50
        
        self.checkBox_handlaunched.clicked.connect(self.refresh_window)
        
        self.pushButton_reject.clicked.connect(self.cancel)
        
        self.doubleSpinBox_payload.valueChanged.connect(lambda: self.payload_changed(shape,polar))
        self.doubleSpinBox_theta_sweep.valueChanged.connect(lambda: self.theta_sweep_changed(shape,polar))
        self.doubleSpinBox_body_length.valueChanged.connect(lambda: self.body_length_changed(shape,polar))
        self.doubleSpinBox_body_height.valueChanged.connect(lambda: self.body_height_changed(shape,polar))
        
        self.verticalSlider_structure_factor.valueChanged.connect(lambda: self.struct_slider_changed(shape,polar))
        self.verticalSlider_aero_factor.valueChanged.connect(lambda: self.aero_slider_changed(shape,polar))
        self.horizontalSlider_span.valueChanged.connect(lambda: self.span_slider_changed(shape,polar))
        self.horizontalSlider_Q_factor.valueChanged.connect(lambda: self.Q_slider_changed(shape,polar))
#        self.doubleSpinBox_body_length.valueChanged.connect(lambda: self.plot_polar(polar))

    def cancel(self):
        self.close()

    def payload_changed(self, shape, polar):
        
        total_mass=self.doubleSpinBox_payload.value()/1000/self.payload_ratio
        self.label_total_mass.setText('{:0.2f} kg'.format(total_mass))
        self.total_mass=total_mass
        
        self.recalc_mass_S()
        self.recalc_span()
        self.draw_plane(shape)
        self.plot_polar(polar)
        
    def body_length_changed(self, shape, polar):
        self.c_root=self.doubleSpinBox_body_length.value()/1000
        
#        b_ratio=self.horizontalSlider_span.value()/100
#        b=self.b_min + b_ratio*(self.b_max-self.b_min)
#        self.b=b
#        self.AR=b**2/self.S_ref
#        self.recalc_mass_S()
        self.recalc_span()
        self.draw_plane(shape)
        self.plot_polar(polar)
        
    def body_height_changed(self, shape, polar):
        self.body_height=self.doubleSpinBox_body_height.value()/1000
        
#        self.recalc_mass_S()
#        self.recalc_span()
#        self.draw_plane(shape)
        self.plot_polar(polar)
        
    def theta_sweep_changed(self, shape, polar):
#        self.recalc_mass_S()
        self.theta_sweep=self.doubleSpinBox_theta_sweep.value()*np.pi/180
        
        self.recalc_span()
        self.draw_plane(shape)
        self.plot_polar(polar)
        
    def struct_slider_changed(self, shape, polar):
        structure_factor=self.verticalSlider_structure_factor.value()/50
        self.payload_ratio=self.payload_ratio_min + structure_factor*(self.payload_ratio_max-self.payload_ratio_min)        
        
        total_mass=self.doubleSpinBox_payload.value()/1000/self.payload_ratio
        self.label_total_mass.setText('{:0.2f} kg'.format(total_mass))
        self.total_mass=total_mass
        
        self.recalc_mass_S()
        self.recalc_span()
        self.draw_plane(shape)
        self.plot_polar(polar)
        
    def aero_slider_changed(self, shape, polar):
        self.aero_factor=self.verticalSlider_aero_factor.value()/50
        self.tc_ratio=0.3 - self.aero_factor*0.2
        self.recalc_mass_S()
        self.recalc_span()
        self.draw_plane(shape)
        self.plot_polar(polar)
    
    def span_slider_changed(self, shape, polar):
        self.b_ratio=self.horizontalSlider_span.value()/100

        self.recalc_span()
        self.draw_plane(shape)
        self.plot_polar(polar)
        
    def Q_slider_changed(self, shape, polar):
        Q_ratio=self.horizontalSlider_Q_factor.value()/50
        self.Q=1.1 + Q_ratio*0.2 #ranges from 1.1 to 1.3
#        self.recalc_mass_S()
        self.recalc_span()
        self.draw_plane(shape)
        self.plot_polar(polar)
        
    def refresh_window(self):
        if self.checkBox_handlaunched.isChecked():
            self.doubleSpinBox_stall_speed.setEnabled(False)

        else:
            self.doubleSpinBox_stall_speed.setEnabled(True)
    
    def recalc_mass_S(self):
        if self.doubleSpinBox_payload.value() == 0:
            return
        try:
            total_mass=self.total_mass
            theta_sweep=self.theta_sweep
            rho=self.rho
        except AttributeError:
            return
        
        CL_max_max=1.4
        CL_max_min=0.6
        
        L=9.8*total_mass
        CL_max=CL_max_min+self.aero_factor*(CL_max_max-CL_max_min) #range from 1 to 2 (best possible with flaps)

        CL_max_swept=CL_max*np.cos(theta_sweep)
        
#        U_thrown_min=4 #Paper on development of portable drones estimates adult male can launch drone at 8m/s
#        static_acc_min=8
#        allowable_drop=0.8
#        
#        drop_time=np.sqrt(2*allowable_drop/9.81)
        U_stall=8#U_thrown_min+static_acc_min*drop_time 8 is about the minimum for one with a decent enough motor, the drop isn't at 9.8m/s anyway and people launch at all weird angles
        
#        print(U_stall)
        
        S_ref=2*L/(rho*U_stall**2*CL_max_swept)
        self.label_wing_area.setText('{:0.2f} m2'.format(S_ref))
        self.S_ref=S_ref
        
        
        
        
    def recalc_span(self):
        
        try:
            S_ref=self.S_ref
            c_root=self.c_root
            b_ratio=self.b_ratio
        except AttributeError:
            return
        
        
        AR_max=4*S_ref/c_root**2 #50.0
    
        AR_min=S_ref/c_root**2 #2
        b_min=np.sqrt(S_ref*AR_min)
        b_max=np.sqrt(S_ref*AR_max)
        
        
        
        b=b_min + b_ratio*(b_max-b_min)
        self.b=b
        
        c_tip=2*S_ref/b-c_root
        self.c_tip=c_tip
        
        self.label_b_min.setText('{:0.2f} m'.format(b_min))
        self.label_b_max.setText('{:0.2f} m'.format(b_max))
        self.label_b_chosen.setText('{:0.2f} m'.format(b))
        self.label_span.setText('{:0.2f} m'.format(b))
        
#        b=self.b_min + self.b_ratio*(self.b_max-self.b_min)
#        self.b=b
        AR=b**2/self.S_ref
        self.AR=AR
        self.label_AR.setText('{:0.2f} m'.format(AR))
#        self.draw_plane(axes)
    
    
    def plot_polar(self, axes):
#        axes.clear()
#        print(self.polar_line)
        try:
            theta_sweep=self.theta_sweep
            rho=self.rho
            S_ref=self.S_ref
            c_root=self.c_root
            AR=self.AR
            b=self.b
            total_mass=self.total_mass
            U=self.Urange
            Q=self.Q
            body_height=self.body_height
            tc_ratio=self.tc_ratio
            
            c_tip=self.c_tip
        except AttributeError:
            return
        
        
#        tc_ratio=0.2
        
        Z=2*np.cos(theta_sweep)
        K=1+Z*(tc_ratio)+100*(tc_ratio)**4 #Form factor of wing, Shevell 1989
        
#        print(K)
        Re=U*np.sqrt(S_ref)*rho/self.atmosphere['mu'] #Reynolds number
        C_f=0.455/(np.log10(Re)**2.58) #Schlichting Prandtl skin friction    
        S_wet=S_ref*2*1.02    
        
        Cd_p=Q/S_ref*K*C_f*S_wet
#        print(Cd_p)
        #Add for fuselage and vertical surfaces
        
        lambda_t=c_tip/c_root - 0.357 + 0.45*np.exp(0.0375*theta_sweep) #Taper ratio with sweep correction
        func_lambda=0.0524*lambda_t**4 - 0.15*lambda_t**3 + 0.1659*lambda_t**2 - 0.0706*lambda_t + 0.0119
        e_th=1/(1+func_lambda*AR)
        
        #Correction factors for the Oswald form factor
        ke_D0=0.9 #Rough estimation for general aviation
        ke_WL=1.00 #Should be accounted for properly maybe
        
        ke_f= 1-2*(body_height/b) #Accounts for fuselage
        
        e_total=e_th*ke_f*ke_D0*ke_WL #total Oswald form factor
#        print(e_total)
        L=9.8*total_mass
        Cl=2*L/(rho*U**2*S_ref)
        
        Cd_if=1/(np.pi*e_total*AR)

        C_d=Cd_p + Cd_if*Cl**2
    #    C_d=0.0174 + 0.1156*Cl**2
    
        self.Cl=Cl
        self.Cd=C_d
    #    Re=U*np.sqrt(plane['ref_area'])*atmosphere['rho']/atmosphere['mu'] #Reynolds number
        drag=C_d*0.5*rho*S_ref*U**2
        drag_min_vel=U[np.argmin(drag)]
        self.label_drag_min_vel.setText('{:0.2f} m/s'.format(drag_min_vel))

#        axes.axis([0, 30, 0, 5*drag_min])        
        
        self.polar_line.set_ydata(drag)
#        axes.draw_artist(axes.yaxis)
        axes.draw_artist(axes.patch)
        axes.draw_artist(self.polar_line)
        self.canvas_polar.update()
        self.canvas_polar.flush_events()
    
#        self.canvas_polar.draw()
        
    def draw_plane(self, axes):
#        axes.clear()

        try:
            b=self.b
            c_root=self.c_root
            c_tip=self.c_tip
            theta_sweep=self.theta_sweep
        except AttributeError:
            return
        spacer=c_root/6
        
        p1y=-0.25*c_root-b/2*np.tan(theta_sweep)+c_tip*0.25
        
        p0x=-b/2
        p0y=p1y-c_tip
        
        p1x=-b/2
        
                
        p2x=0
        p2y=0#c_root+spacer
        
        p3x=b/2
        p3y=p1y

        p4x=b/2
        p4y=p0y  
        
        p5x=0
        p5y=-c_root
        
        verts = [
        (p0x, p0y), # left, bottom
        (p1x, p1y), # left, top
        (p2x, p2y), # right, top
        (p3x, p3y), # right, bottom
        (p4x, p4y), # ignored
        (p5x, p5y), # ignored
        ]
    
#        codes = [Path.MOVETO,
#                 Path.LINETO,
#                 Path.LINETO,
#                 Path.LINETO,
#                 Path.LINETO,
#                 Path.LINETO,
#                 Path.CLOSEPOLY,
#                 ]
        
#        path = Path(verts, codes)
        
#        patch = patches.PathPatch(path, facecolor='gray', lw=2)
#        self.plane_patch=axes.add_patch(patch)
#        print(verts)
#        axes.axis('off')
#        self.plane_patch=patch
        
#        axes.draw_artist(axes.yaxis)
#        axes.axis([0, b/2*1.1, min([0,p3y-spacer]), p1y+spacer])
        if abs(p3x-self.xlim_old)>spacer:
#            axes.axis('equal')
#            self.canvas_shape.restore_region(self.background)
            axes.set_xlim([-p3x-spacer,p3x+spacer])
            distance=2*(p3x+spacer)
            centre=-spacer-c_root/2
            axes.set_ylim([centre-distance/2,centre+distance/2])
            self.xlim_old=p3x
#            axes.draw_artist(axes.yaxis)
#            axes.draw_artist(axes.xaxis)
            self.canvas_shape.draw()
            self.background = self.canvas_shape.copy_from_bbox(axes.bbox)
#        axes.draw_artist(axes.patch)
#        self.background = self.canvas_shape.copy_from_bbox(axes.bbox)
        self.plane_patch.set_xy(verts)
        self.canvas_shape.restore_region(self.background)
        axes.draw_artist(self.plane_patch)
#        self.canvas_shape.blit(axes.bbox)
#        self.background = self.canvas_shape.copy_from_bbox(axes.bbox)
        
#        axes.draw_artist(axes.patch)
        self.canvas_shape.update()
        self.canvas_shape.flush_events()
        
#        self.canvas_shape.draw()