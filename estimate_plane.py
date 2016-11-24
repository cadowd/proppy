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


class estimate_plane_window(QtGui.QDialog):
    def __init__(self):
        super(estimate_plane_window, self).__init__()
        uic.loadUi("./UI/estimator_dialog.ui", self)
        self.pushButton_calc_props.clicked.connect(self.calc_props)
    
    def refresh(self):
        if self.checkBox_hand_launched.isChecked():
            self.doubleSpinBox_stall_speed.setEnabled(False)

        else:
            self.doubleSpinBox_stall_speed.setEnabled(True)
