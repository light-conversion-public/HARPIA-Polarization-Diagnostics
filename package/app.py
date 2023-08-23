#!/usr/bin/env python
# -*- coding: utf-8 -*-
#==========================================================================
# HARPIA Polarization Diagnostics Tool
#--------------------------------------------------------------------------
# Copyright (c) 2022 Light Conversion, UAB
# All rights reserved.
# www.lightcon.com
#==========================================================================
     
import lclauncher

import os
import sys    
import time
import numpy as np
import lightcon.style
import json

# lightcon.style.apply_style()

sys.path.append(os.path.dirname(os.path.realpath(sys.argv[0])))
os.chdir(os.path.dirname(os.path.realpath(sys.argv[0])))

# if connections to devices are used, they are initiated here:
connections = lclauncher.establish_connections()

# initialize and connect to HARPIA
harpia = connections.get_connection('harpia')

# check if connection successful
if not harpia:
    sys.exit("Could not connect to Harpia")

import matplotlib
matplotlib.use('Qt5Agg')

from PyQt5.QtWidgets import *
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot, QObject, QThread, pyqtSignal
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure


from utils import *

f = open('./package/settings.json', 'r')
settings = json.loads(f.read())

motor_index = settings["motor_index"]
can_id = settings["can_id"]
reduction = settings["reduction"]

old_spectra_per_acquisition = harpia.spectra_per_acquisition()

can_sender = HarpiaCanSender(harpia)

lcan = LepreCanDevice(None, can_id)
mb = MotorBoard(can_id, can_sender)

diag_unit = PolarizationDiagnosticsUnit(mb, motor_index, reduction, settings.get("speed") or 10000, r'package/Sanyo Denki SH2281-5631 (rotary).json', settings.get("zero_angle") or 0.0)

angles = np.array([])
intensities = np.array([])

def get_intensity():
    return np.abs(np.average(harpia._get('Basic/RawSignal')[settings['signal']]))

def analyze(intensities):
    if len(intensities) == 0:
        return None
    
    t1_ind = np.argmax(intensities)
    T1 = np.max(intensities)
    t2_ind = np.argmin(intensities)
    T2 = intensities[t2_ind]
    
    P = (T1 - T2) / (T1 + T2)
    rho_P = T2/T1
    
    return {"extinction": 1.0/rho_P, "efficiency": P, "min_ind": t2_ind}

class Worker(QObject):
    finished = pyqtSignal()
    progress = pyqtSignal(dict)

    def run(self):
        self.is_running = True
        
        angles = np.array([])
        intensities = np.array([])
        
        harpia.set_spectra_per_acquisition(settings['spectra_per_acquisition'])
                
        diag_unit.stop()
        
        init_angle = diag_unit.get_angle() % 360.0
        current_angle = init_angle
        
        diag_unit.start()   

        cnt = 0.0
        while cnt < 360.0:
            old_angle = current_angle
            current_angle = diag_unit.get_angle()
            if current_angle < old_angle:
                cnt = cnt + current_angle + 360.0 - old_angle
            else:
                cnt = cnt + current_angle - old_angle
            angles = np.append(angles, [(current_angle) / 180.0 * np.pi])
            intensities = np.append(intensities, [get_intensity()])
            self.progress.emit({'finalize': False, 'angles': angles, 'intensities': intensities})
        
        self.progress.emit({'finalize': True, 'angles': angles, 'intensities': intensities})
                
        diag_unit.stop()
        
        diag_unit.set_angle_blocking(0.0)
        
        harpia.set_spectra_per_acquisition(old_spectra_per_acquisition)

        self.finished.emit()
        
    def reset(self):
        self.is_running = True
        
        diag_unit.reset()        
    
        self.finished.emit()
        
class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        lightcon.style.apply_style()
        self.fig = Figure(figsize=(width, height), dpi=dpi)                
        
        plt.ion()
        # plt.clf()        
        
        self.ax = self.fig.add_subplot(111, projection = 'polar')

        self.fig.suptitle("Transmission")

        self.line, = self.ax.plot([], [], label='actual')
        self.maxpolline, = self.ax.plot([0, np.pi], [10.0] * 2, 'k--', lw = 1.0)
        self.minpolline, = self.ax.plot([np.pi/2.0, np.pi + np.pi/2.0], [10.0] * 2, 'k--', lw = 1.0)

        self.ax.set_ylim([0.0, 10.0])
        self.ax.set_thetalim(0, 2.0*np.pi)
        self.ax.set_xticks(np.arange(0.0, 360.0, 10) / 180.0 * np.pi)
        self.ax.yaxis.set_ticklabels([])
        
        lightcon.style.add_watermark(self.ax)
                
        super(MplCanvas, self).__init__(self.fig)
        
class MainWindow(QMainWindow):
    sc = None
    canDraw = True
    
    def __init__(self, title):
        super().__init__()
        self.title = title
        self.left = 100
        self.top = 100
        self.width = 900
        self.height = 900
        self.initUI()
        
            
    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        
        self.start_button = QPushButton('START', self)
        self.start_button.setToolTip('Start measuring')
        self.start_button.clicked.connect(self.start_button_on_click)

        self.reset_button = QPushButton('RESET', self)
        self.reset_button.setToolTip('Reset unit')
        self.reset_button.clicked.connect(self.reset_button_on_click)

        # self.clear_button = QPushButton('CLEAR', self)
        # self.clear_button.clicked.connect(self.clear_button_on_click)
                
        self.sc = [MplCanvas(self, width=5, height=3, dpi=100)]

        outerLayout = QVBoxLayout()
        topLayout = QHBoxLayout()            
        bottomLayout = QVBoxLayout()

        topLayout.addWidget(self.start_button, 1)
        topLayout.addWidget(self.reset_button, 1)
        # topLayout.addWidget(self.clear_button, 1)

        bottomLayout.addWidget(self.sc[0])        
        
        outerLayout.addLayout(topLayout)
        outerLayout.addLayout(bottomLayout)
        
        widget = QWidget()
        widget.setLayout(outerLayout)

        self.setCentralWidget(widget)  
        self.show()
        
    def runLongTask(self):        
        self.thread = QThread()
        self.worker = Worker()
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)
        self.thread.finished.connect(self.enableButtons)
        self.worker.progress.connect(self.updatePlots)
        self.thread.start()
        # Final resets       
        # self.thread.finished.connect(self.addToPlots)

    def resetLongTask(self):
        self.thread = QThread()
        self.worker = Worker()
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.reset)
        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)
        self.thread.finished.connect(self.enableButtons)        
        self.thread.start()

    def updatePlots(self, info): 
        if (not self.canDraw) and (not info['finalize']): return()
        self.canDraw = False
        intensities = info['intensities']
        angles = info['angles']
        if info['finalize'] :
            self.sc[0].line.set_xdata(np.concatenate(((angles[angles<np.pi*2])[:-1], [angles[0] + np.pi * 2])))
            self.sc[0].line.set_ydata(np.concatenate(((intensities[angles<np.pi*2])[:-1], [intensities[0]])))
        else:            
            self.sc[0].line.set_xdata((angles[angles<np.pi*2])[:-1])
            self.sc[0].line.set_ydata((intensities[angles<np.pi*2])[:-1])
        
            pars = analyze(intensities)        
            if pars:
                tangle = angles[pars["min_ind"]]
                self.sc[0].maxpolline.set_xdata(np.array([tangle, np.pi+tangle]))
                self.sc[0].maxpolline.set_ydata([np.max(intensities)]*2)
                self.sc[0].minpolline.set_xdata([tangle + np.pi/2.0, np.pi + np.pi/2.0 +tangle])
                self.sc[0].minpolline.set_ydata([np.max(intensities)]*2)
            
                self.sc[0].ax.set_title("{:.2f}:1, efficiency {:.2f}".format(pars["extinction"], pars["efficiency"]))
        
                self.sc[0].ax.set_ylim([0, np.max(intensities)])
        
        self.sc[0].draw()
        app.processEvents()
        self.canDraw = True
        
        
        # if (len(self.beam_log) > 0):
        #     if info['position'] == self.beam_log[-1]['position']:
        #         self.beam_log[-1] = info
        #     else:
        #         self.beam_log = self.beam_log + [info]
        # else:
        #     self.beam_log = self.beam_log + [info]

        # imin = np.argmin([item['position'] for item in self.beam_log])
        # imax = np.argmax([item['position'] for item in self.beam_log])

        # if (imin < imax):
        #     if (imax != len(self.beam_log) - 1):
        #         self.beam_log = self.beam_log[imax:]
        # else:
        #     if (imin != len(self.beam_log) - 1):
        #         self.beam_log = self.beam_log[imin:]
        
        # # self.beam_log = self.beam_log[-100:]

        # for sc in self.sc:
        #     sc.ax1.cla()
        #     sc.ax2.cla()
        #     sc.ax1.set_xlabel('Delay line position, mm')

        # x = [item['position'] for item in self.beam_log]
        # fx = [x[0],x[-1]]
        # y = [[item['beam_parameters']['MeanX'] for item in self.beam_log], [item['beam_parameters']['MeanY'] for item in self.beam_log]]
        # f = [np.poly1d(np.polyfit(x,y[0],1))(fx), np.poly1d(np.polyfit(x,y[1],1))(fx) ]
        # lns = []

        # lns = lns + self.sc[0].ax1.plot(x, y[0], '.-', color='C0', label = '$\\Delta X$ {:.0f} um'.format(np.max(y[0]) - np.min(y[0])))
        # lns = lns + self.sc[0].ax2.plot(x, y[1], '.-',color='C1', label = '$\\Delta Y$ {:.0f} um'.format(np.max(y[1]) - np.min(y[1])))    
        # self.sc[0].ax1.plot(fx, f[0], '--', color='C0')
        # self.sc[0].ax2.plot(fx, f[1], '--',color='C1')  
        # self.sc[0].ax1.set_ylabel('X beam position, mm')
        # self.sc[0].ax2.set_ylabel('Y beam position, mm')
        
        # self.sc[0].ax1.legend(lns, [l.get_label() for l in lns])

        # y = [[item['beam_parameters']['SigmaPrimary'] for item in self.beam_log], [item['beam_parameters']['SigmaSecondary'] for item in self.beam_log]]
        # f = [np.poly1d(np.polyfit(x,y[0],1))(fx), np.poly1d(np.polyfit(x,y[1],1))(fx) ]

        # self.sc[1].ax1.plot(x, y[0], '.-', color='C0')
        # self.sc[1].ax2.plot(x, y[1], '.-',color='C1')
        # self.sc[1].ax1.plot(fx, f[0], '--', color='C0')
        # self.sc[1].ax2.plot(fx, f[1], '--',color='C1')  
        # self.sc[1].ax1.set_ylabel('X beam sigma P, mm')
        # self.sc[1].ax2.set_ylabel('Y beam sigma S, mm')

        # for sc in self.sc:
        #     sc.draw_idle()
    def enableButtons(self):
        self.reset_button.setEnabled(True)
        self.start_button.setEnabled(True)
    def disableButtons(self):
        self.reset_button.setEnabled(False)
        self.start_button.setEnabled(False)
        
    @pyqtSlot()
    def start_button_on_click(self):
        # camera.enable_beam_profiler()
        self.runLongTask()
        self.disableButtons()
        
    @pyqtSlot()
    def reset_button_on_click(self):
        # camera.enable_beam_profiler()
        self.resetLongTask()
        self.disableButtons()


    @pyqtSlot()
    def clear_button_on_click(self):
        # camera.enable_beam_profiler()
        self.beam_log = []

app = QApplication([])
w = MainWindow('HARPIA Polarization diagnostics')
app.exec_()
