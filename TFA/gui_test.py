# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 09:09:54 2017

@author: NTU_Math
"""

import sys
import numpy as np
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
from functions import *
from CFPH_functions import *
from wvd_fun import *
import matplotlib.colors as colors
from matplotlib.widgets import Slider

x_raw = []; t = []; t_todo = np.arange(0,100); signal_todo = []
plot_count = 2
total_plot = 6
if_plot = [False, True, True, False, False, False]
if_log = [0,0,0,0,0,0]
plot_row_n = 1
SR = 100 
view_wl = 100
plot0_data = []
plot1_data = []
plot2_data = []
plot3_data = []
plot4_data = []
plot5_data = []


class Window(QMainWindow):
    
    
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)
        """Menu Bar"""
        bar = self.menuBar()
        fileMenu = bar.addMenu('&File')
        fileMenu.addAction("Open Data")
        fileMenu.triggered[QAction].connect(self.windowaction)        
        plotMenu = bar.addMenu('Plot')
        plotMenu.addAction("Plot Selection")
        plotMenu.triggered[QAction].connect(self.windowaction)
        self.form_widget = FormWidget(self) 
        self.form_widget.run_plot()
        self.setCentralWidget(self.form_widget)
        

        """Window Title"""
        self.setWindowTitle("Time Frequency Analysis")
		
    def windowaction(self, q):
        
        if q.text() == "Open Data":
            global x_raw, t, t_todo, signal_todo, plot1_data, plot2_data, plot3_data, plot4_data, plot5_data
            plot0_data = []
            plot1_data = []
            plot2_data = []
            plot3_data = []
            plot4_data = []
            plot5_data = []
            self.name = QFileDialog.getOpenFileName(self, 'Open Data','' ,"*.txt *.csv")
            
            if self.name[0]:
                data = np.genfromtxt(self.name[0], delimiter=',')
                n = len(data.shape)
                if n == 1:
                    x_raw = data; t = np.arange(0,len(data))
                    t_todo = t; signal_todo = x_raw
                    
                elif n==2:
                    t = data[:,0]; x_raw = data[:,1]
                    t_todo = t; signal_todo = x_raw
                else:
                    print ('Wrong Format')
        
        elif q.text() == "Plot Selection":
            self.ps_dialog = QDialog()
            self.plot_cb = []
            for i in range(0,total_plot):
                self.plot_cb.append(QCheckBox())
            
            self.plot_cb[0].setText('Power Spectrum')
            self.plot_cb[1].setText('Signal View')
            self.plot_cb[2].setText('STFT')
            self.plot_cb[3].setText('CWT')
            self.plot_cb[4].setText('S Transform')
            self.plot_cb[5].setText('WVD')
            
            
            
            vbox = QFormLayout(self.ps_dialog)
            self.plot_cb[0].setChecked(if_plot[0])
            vbox.addRow('Plot Selection', self.plot_cb[0])
            for i in range(1,total_plot):
                self.plot_cb[i].setChecked(if_plot[i])
                vbox.addWidget(self.plot_cb[i])
            
            self.PlotRow = QLineEdit(str(plot_row_n))
            self.PlotRow.setValidator(QIntValidator())
            

            Done = QPushButton('Done')
            Done.clicked.connect(self.plot_selection_Done)
            vbox.addRow('Number of Rows in plot', self.PlotRow)
            vbox.addWidget(Done)
            self.setLayout(vbox)
            self.ps_dialog.setWindowTitle("Plot Selection")
            self.ps_dialog.setWindowModality(Qt.ApplicationModal)
            self.ps_dialog.exec_()
            
    
    def plot_selection_Done(self):
        global if_plot, plot_count, plot_row_n
        
        for i in range(0,total_plot):
            if_plot[i] = self.plot_cb[i].isChecked()
            
        plot_row_n = int(self.PlotRow.text())
        plot_count = np.sum(if_plot)
        self.form_widget.run_plot()
        self.ps_dialog.accept()
        
        
        
        
    
class FormWidget(QWidget):

    def __init__(self, parent):        
        super(FormWidget, self).__init__(parent)
        #Figure1
        
        self.figure_hide = plt.figure()
        self.ax_hide = self.figure_hide.add_subplot('111')
        self.mw_layout_vbox = QVBoxLayout()
        
        e_t = np.arange(0,100)
        e_y = np.arange(0,100)
        e_M = np.zeros((100,100))
        e_M[0,0] = 10           
        self.e_cax = self.ax_hide.pcolorfast(e_t, e_y, e_M, cmap = 'Greys')
        self.ax_hide.set_xlim([0,0+view_wl])
        
        #slider
        self.sl = QSlider(Qt.Horizontal)
        self.sl_min = QLabel(str(np.min(t_todo)))
        self.sl.setMinimum(np.min(t_todo))
        self.sl_max = QLabel(str(np.max(t_todo)-view_wl))
        self.sl.setMaximum(np.max(t_todo)-view_wl)
        self.sl.setValue(np.min(t_todo))
        self.sl_cur_val = QLabel(str(self.sl.value()))
        self.sl.setTickPosition(QSlider.TicksBelow)
        self.sl.setTickInterval(5)
        self.sl.valueChanged.connect(self.sl_valuechange)
        self.sl_wl = QLineEdit(str(view_wl))
        self.sl_wl.setFixedWidth(50)
        self.sl_wl.setValidator(QIntValidator())
        
        
        
        


             
    def clearLayout(self, layout):
        for i in reversed(range(layout.count())):
            item = layout.itemAt(i)

            if isinstance(item, QWidgetItem):
                print ("widget" + str(item))
                item.widget().close()
                # or
                # item.widget().setParent(None)
            elif isinstance(item, QSpacerItem):
                print ("spacer " + str(item))
                # no need to do extra stuff
            else:
                print ("layout " + str(item))
                self.clearLayout(item.layout())

            # remove the item from layout
            layout.removeItem(item)    
    
    def run_plot(self):
        self.clearLayout(self.mw_layout_vbox)
        
        order = 0
        self.mw_layout_hbox = QHBoxLayout()
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.panel_base = QWidget()
        self.panel = QVBoxLayout()
        if plot_count>3 and plot_row_n<2:
            self.nrow = 2
        else:
            self.nrow = plot_row_n
        self.ncol = int(np.ceil(plot_count/self.nrow))
        for i in reversed(range(self.panel.count())): 
            self.panel.itemAt(i).widget().deleteLater()
            
        self.gb_pre = QGroupBox('Preprocessing')
        self.para_layout_pre = QVBoxLayout()
        self.filter_btype = QComboBox()
        self.filter_btype.addItems(['None','lowpass','highpass','bandpass'])
        self.button_pre = QPushButton('Update')
        self.button_pre.clicked.connect(self.update_pre)
        self.para_layout_pre.addWidget(self.filter_btype)
        self.para_layout_pre.addWidget(self.button_pre)
        self.gb_pre.setLayout(self.para_layout_pre)
        self.panel.addWidget(self.gb_pre)
        
        
        if if_plot[0]:
            #Figure0
            order+=1
            self.ax0 = self.figure.add_subplot(str(self.nrow)+str(self.ncol)+str(order))
            self.ax0.grid(linestyle='dotted')
            if plot0_data:
                self.ax0.semilogy(plot0_data[0], plot0_data[1],color = 'k',linewidth = 0.3)
            self.ax0.set_title('Power Spectrum')
            self.ax0.set_xlabel('frequency [Hz]')
            self.ax0.set_ylabel('V**2')
            self.gb_ps = QGroupBox('Power Spectrum')
            self.para_layout0 = QVBoxLayout()
            self.button_0 = QPushButton('Update')
            self.button_0.clicked.connect(lambda:self.update_ps(self.ax0))
            self.para_layout0.addWidget(self.button_0)
            self.gb_ps.setLayout(self.para_layout0)
            self.panel.addWidget(self.gb_ps)
            
                
            
        if if_plot[1]:
            #Figure1
            order+=1
            
            self.ax1 = self.figure.add_subplot(str(self.nrow)+str(self.ncol)+str(order),sharex = self.ax_hide)
            self.ax1.grid(linestyle='dotted')
            if plot1_data:
                self.ax1.plot(plot1_data[0],plot1_data[1], '-', color = 'k', linewidth = 0.3)
                self.ax1.set_xlim([np.min(plot1_data[0]),np.max(plot1_data[0])])
            self.ax1.set_title('Signal View')
            self.ax1.set_xlabel('Time')
            self.gb_signal = QGroupBox('Signal View')
            self.para_layout1 = QVBoxLayout()
            self.button_1 = QPushButton('Update')
            self.button_1.clicked.connect(self.update_1)
            self.para_layout1.addWidget(self.button_1)
            self.gb_signal.setLayout(self.para_layout1)
            self.panel.addWidget(self.gb_signal)
                
        
        if if_plot[2]:
            #Figure2
            order+=1
            self.ax2 = self.figure.add_subplot(str(self.nrow)+str(self.ncol)+str(order),sharex = self.ax_hide,sharey=self.ax_hide)
            self.ax2.grid(linestyle='dotted')
            
            if not plot2_data:
                self.colorbar_2 = self.figure.colorbar(self.e_cax, ax=self.ax2, extend='max')
            else:
                self.colorbar_2 = imageSQ(self.figure, self.ax2, self.colorbar_2, plot2_data[0], plot2_data[1]*SR, np.abs(plot2_data[2]), 99.5,if_log[1], 'Greys')
            self.ax2.set_title('STFT')
            self.ax2.set_xlabel('Time')
            self.ax2.set_ylabel('Frequency')
            self.ax2_signal = self.ax2.twinx()

            self.gb_STFT = QGroupBox('STFT')
            self.para_layout2 = QVBoxLayout()
            self.zlog_cb2 = QComboBox()
            self.zlog_cb2.addItems(["normal scale","log scale"])
            self.en_cb_2 = QComboBox()
            self.en_cb_2.addItems(["Normal", "Synchrosquezzing", "Reassigment"])
            self.button_2 = QPushButton('Update')
            self.button_2.clicked.connect(self.update_2)
            self.para_layout2.addWidget(self.zlog_cb2)
            self.para_layout2.addWidget(self.en_cb_2)
            self.para_layout2.addWidget(self.button_2)
            self.gb_STFT.setLayout(self.para_layout2)
            self.panel.addWidget(self.gb_STFT)
            
        
        if if_plot[3]:
            #Figure 3
            order+=1
            self.ax3 = self.figure.add_subplot(str(self.nrow)+str(self.ncol)+str(order),sharex = self.ax_hide,sharey=self.ax_hide)
            self.ax3.grid(linestyle='dotted')
            if not plot3_data:
                self.colorbar_3 = self.figure.colorbar(self.e_cax, ax=self.ax3, extend='max')
            else:
                self.colorbar_3 = imageSQ(self.figure, self.ax3, self.colorbar_3, plot3_data[0], plot3_data[1]*SR, np.abs(plot3_data[2]), 99.5,if_log[2], 'Greys')
            self.ax3.set_title('Continous Wavelet Transform')
            self.ax3.set_xlabel('Time')
            self.ax3.set_ylabel('Frequency')
            
            self.gb_CWT = QGroupBox('CWT')
            self.para_layout3 = QVBoxLayout()
            self.zlog_cb3 = QComboBox()
            self.zlog_cb3.addItems(["normal scale","log scale"])
            self.en_cb_3 = QComboBox()
            self.en_cb_3.addItems(["Normal", "Synchrosquezzing"])
            self.MW_cb3 = QComboBox()
            self.MW_cb3.addItems(['Cinfc','morse', 'morse-a', 'morse-b', 'morse-c', 'morlet', 'gaussian', 'meyer', 'BL3'])
            self.button_3 = QPushButton('Update')
            self.button_3.clicked.connect(self.update_3)
            self.para_layout3.addWidget(self.zlog_cb3)
            self.para_layout3.addWidget(self.en_cb_3)
            self.para_layout3.addWidget(self.MW_cb3)
            self.para_layout3.addWidget(self.button_3)
            self.gb_CWT.setLayout(self.para_layout3)
            self.panel.addWidget(self.gb_CWT)
            
        
        if if_plot[4]:
            #Figure 4 
            order+=1
            self.ax4 = self.figure.add_subplot(str(self.nrow)+str(self.ncol)+str(order),sharex = self.ax_hide,sharey=self.ax_hide)
            self.ax4.grid(linestyle='dotted')
            if not plot4_data:
                self.colorbar_4 = self.figure.colorbar(self.e_cax, ax=self.ax4, extend='max')
            else:
                self.colorbar_4 = imageSQ(self.figure, self.ax4, self.colorbar_4, plot4_data[0], plot4_data[1]*SR, np.abs(plot4_data[2]), 99.5,if_log[3], 'Greys')
            self.ax4 = self.figure.add_subplot(str(self.nrow)+str(self.ncol)+str(order),sharex = self.ax_hide,sharey=self.ax_hide)
            self.ax4.set_title('S Transform')
            self.ax4.set_xlabel('Time')
            self.ax4.set_ylabel('Frequency')
            self.gb_stran = QGroupBox('S Transform')
            self.para_layout4 = QVBoxLayout()
            self.zlog_cb4 = QComboBox()
            self.zlog_cb4.addItems(["normal scale","log scale"])
            self.en_cb_4 = QComboBox()
            self.en_cb_4.addItems(["Normal", "Synchrosquezzing"])
            self.button_4 = QPushButton('Update')
            self.button_4.clicked.connect(self.update_4)
            self.para_layout4.addWidget(self.zlog_cb4)
            self.para_layout4.addWidget(self.en_cb_4)
            self.para_layout4.addWidget(self.button_4)
            self.gb_stran.setLayout(self.para_layout4)
            self.panel.addWidget(self.gb_stran)
            
        
        
        if if_plot[5]:
            #Figure5
            order+=1
            self.ax5 = self.figure.add_subplot(str(self.nrow)+str(self.ncol)+str(order),sharex = self.ax_hide,sharey=self.ax_hide)
            self.ax5.grid(linestyle='dotted')
            if not plot5_data:
                self.colorbar_5 = self.figure.colorbar(self.e_cax, ax=self.ax5, extend='max')
            else:
                self.colorbar_5 = imageSQ(self.figure, self.ax5, self.colorbar_5, plot5_data[0], plot5_data[1]*SR, np.abs(plot5_data[2]), 99.5,if_log[4], 'Greys')
            
            self.ax5.set_title('WVD')
            self.ax5.set_xlabel('Time')
            self.ax5.set_ylabel('Frequency')
            self.gb_wvd = QGroupBox('WVD')
            self.para_layout5 = QVBoxLayout()
            self.zlog_cb5 = QComboBox()
            self.zlog_cb5.addItems(["normal scale","log scale"])
            self.en_cb_5 = QComboBox()
            self.en_cb_5.addItems(["Normal", "Synchrosquezzing"])
            self.button_5 = QPushButton('Update')
            self.button_5.clicked.connect(self.update_5)
            self.para_layout5.addWidget(self.zlog_cb5)
            self.para_layout5.addWidget(self.en_cb_5)
            self.para_layout5.addWidget(self.button_5)
            self.gb_wvd.setLayout(self.para_layout5)
            self.panel.addWidget(self.gb_wvd)
            
            
        
        
        self.sl_panel_base = QWidget()
        self.sl_vlayout = QVBoxLayout()
        self.sl_hlayout1 =QHBoxLayout()
        self.sl_hlayout1.addWidget(self.sl_wl)
        self.sl_hlayout1.addWidget(self.sl_cur_val)
        self.sl_vlayout.addLayout(self.sl_hlayout1)
        self.sl_hlayout2 = QHBoxLayout()
        self.sl_hlayout2.addWidget(self.sl_min)
        self.sl_hlayout2.addWidget(self.sl)
        self.sl_hlayout2.addWidget(self.sl_max)
        self.sl_vlayout.addLayout(self.sl_hlayout2)
        self.sl_panel_base.setFixedHeight(50)
        self.sl_panel_base.setLayout(self.sl_vlayout)
        self.panel_base.setLayout(self.panel)
        self.panel_base.setFixedWidth(150)

        

        # set the layout
        self.figure.tight_layout()
        "self.figure.subplotpars.update(left=0.05, bottom=0.05, right=0.95, top=0.95, wspace=0.2, hspace=0.2)"
        
        
        self.mw_layout_vbox.addWidget(self.toolbar)
        self.mw_layout_hbox.addWidget(self.panel_base)
        self.figure_vbox = QVBoxLayout()
        self.figure_vbox.addWidget(self.canvas)
        self.figure_vbox.addWidget(self.sl_panel_base)
        self.mw_layout_hbox.addLayout(self.figure_vbox)
        self.mw_layout_vbox.addLayout(self.mw_layout_hbox)
        self.setLayout(self.mw_layout_vbox)
        self.sl_valuechange()

    def sl_valuechange(self):
        global view_wl
        view_wl = int(self.sl_wl.text())
        t_min = np.min(t_todo)
        t_max = np.max(t_todo)
        if view_wl>=len(t_todo):
            print('Window Length too Long')
            view_wl = int(len(t_todo))
            self.sl.setMaximum(t_max)
            self.sl.setMinimum(t_min)
            self.sl_max.setText(str(t_max))
            self.sl_min.setText(str(t_min))
            self.sl_cur_val.setText(str(t_min))
            self.ax_hide.set_xlim([np.min(t_todo),np.max(t_todo)])
            self.sl.setValue(t_min)
        else:
            pos = self.sl.value()
            if (pos<t_min) or (pos>t_max-view_wl):
                pos = t_min
                
            self.sl.setValue(pos)    
            self.sl_cur_val.setText(str(pos))
            self.ax_hide.set_xlim([pos,pos+view_wl])
            self.sl.setMaximum(t_max-view_wl)
            self.sl.setMinimum(t_min)
            self.sl_max.setText(str(t_max-view_wl))
            self.sl_min.setText(str(t_min))
    

    
    def update_ps(self,ax):
        def plot_0():
            
            from scipy import signal
            global SR, plot0_data
            SR = int(SamplingRate.text())
            WL = int(WindowLength.text())
            scale = Scaling.currentText()
            x, y = signal.periodogram(signal_todo, SR,nfft = WL, scaling = scale)
            plot0_data = []
            plot0_data.append(x), plot0_data.append(y)
            ax.hold(False)
            ax.semilogy(x[1:], y[1:],  color = 'k',linewidth = 0.3)
            ax.set_title('Power Spectrum')
            ax.set_xlabel('frequency [Hz]')
            ax.set_ylabel('V**2')
            ax.grid(linestyle='dotted')
            self.dialog.accept()
            self.canvas.draw()
        
        if signal_todo == []:
            print ('data not loaded')
            return
        
        
        self.dialog = QDialog()
        vbox = QVBoxLayout(self.dialog)
        #Group Box1 for General Parameeters
        group_box0 = QGroupBox()
        gb0 = QFormLayout()
        
        SamplingRate = QLineEdit(str(SR))
        SamplingRate.setValidator(QIntValidator())
        WindowLength = QLineEdit(str(len(t_todo)))
        WindowLength.setValidator(QIntValidator())
        Scaling = QComboBox()
        Scaling.addItems(['spectrum','density'])
        gb0.addRow("Sampling Rate",SamplingRate)
        gb0.addRow('Window Length', WindowLength)
        gb0.addRow('Scaling', Scaling)
        
        
        group_box0.setLayout(gb0)
        
        
        #Run Button
        Run = QPushButton('Run')
        Run.clicked.connect(plot_0)
        
        vbox.addWidget(group_box0)
        vbox.addWidget(Run)
        
        self.setLayout(vbox)
        self.dialog.setWindowTitle("Power Spectrun Settings")
        self.dialog.setWindowModality(Qt.ApplicationModal)
        self.dialog.exec_()
    
    def update_pre(self):
        def plot_pre():
            global t_todo, signal_todo, SR
            DT = detrend.isChecked()
            SR = int(SamplingRate.text())
            DS = int(DownSample.text())
            N= int(order.text())
            cutoff = []
            if (btype == 'lowpass')or(btype == 'highpass') :
                cutoff = float(cutoff_freq.text())
            elif btype == 'bandpass':
                cutoff = [float(low_freq.text()),float(high_freq.text())]
                if cutoff[0]>=cutoff[1]:
                    print ('Low Freq cannot be higher than High Freq')
                    self.dialog.accept()
                    return
            rp = float(maximum_ripple.text())
            rs = float(minimum_attenuation.text())
            if btype!='None':
                t_todo, signal_todo = signal_filter(t,x_raw,SR,DT,DS, btype,ftype.currentText(),N,cutoff,rp,rs)
            else:
                t_todo, signal_todo = signal_filter(t,x_raw,SR,DT,DS, btype)
            
            self.dialog.accept()
            self.canvas.draw()
        
        
        if signal_todo == []:
            print ('data not loaded')
            return
        
        
        self.dialog = QDialog()
        vbox = QVBoxLayout(self.dialog)
        btype = self.filter_btype.currentText()
        #Group Box1 for General Parameeters
        group_box_pre = QGroupBox()
        gb_pre = QFormLayout()
        
        detrend = QCheckBox('Detrend')
        detrend.setChecked(False)
        SamplingRate = QLineEdit(str(SR))
        SamplingRate.setValidator(QIntValidator())
        DownSample = QLineEdit('0')
        DownSample.setValidator(QIntValidator())
        ftype = QComboBox()
        ftype.addItems(['butter','cheby1','cheby2','ellip','bessel'])
        order = QLineEdit('3')
        order.setValidator(QIntValidator(0,13))
        cutoff_freq = QLineEdit('0.3')
        cutoff_freq.setValidator(QDoubleValidator(0,0.5,3))
        low_freq = QLineEdit('0.1')
        low_freq.setValidator(QDoubleValidator(0,0.5,3))
        high_freq = QLineEdit('0.4')
        high_freq.setValidator(QDoubleValidator(0,0.5,3))
        maximum_ripple = QLineEdit('50')
        maximum_ripple.setValidator(QDoubleValidator())
        minimum_attenuation = QLineEdit('60')
        minimum_attenuation.setValidator(QDoubleValidator())
        
        gb_pre.addWidget(detrend)
        gb_pre.addRow("Sampling Rate",SamplingRate)
        gb_pre.addRow('Down Sample', DownSample)
        
        if btype!='None':
            gb_pre.addRow('Filter Type', ftype)
            gb_pre.addRow('Filter Order', order)
            if (btype == 'lowpass')or(btype == 'highpass') :
                gb_pre.addRow("Cutoff Frequency", cutoff_freq)
            elif btype == 'bandpass':
                gb_pre.addRow("Low Frequency", low_freq)
                gb_pre.addRow("High Frequency", high_freq)
            gb_pre.addRow('Maximum Ripple(dB)', maximum_ripple)
            gb_pre.addRow('Minimum Attenuation(dB)', minimum_attenuation)
        group_box_pre.setLayout(gb_pre)
        
        
        #Run Button
        Run = QPushButton('Run')
        Run.clicked.connect(plot_pre)
        
        vbox.addWidget(group_box_pre)
        vbox.addWidget(Run)
        
        self.setLayout(vbox)
        self.dialog.setWindowTitle("Preprocess Settings")
        self.dialog.setWindowModality(Qt.ApplicationModal)
        self.dialog.exec_()
    
    def update_1(self):
        def plot_1():
            global SR, plot1_data
            SR = int(SamplingRate.text())
            
            plot1_data = []
            plot1_data.append(t_todo), plot1_data.append(signal_todo)
            self.ax1.hold(False)
            self.ax1.plot(t_todo,signal_todo, '-', color = 'k', linewidth = 0.3)
            self.ax1.set_title('Signal View')
            self.ax1.set_xlabel('Time')
            self.ax1.set_xlim([np.min(t_todo),np.max(t_todo)])
            self.ax1.grid(linestyle='dotted')
            self.sl_valuechange()
            self.dialog.accept()
            self.canvas.draw()
        
        if signal_todo == []:
            print ('data not loaded')
            return
        
        
        self.dialog = QDialog()
        vbox = QVBoxLayout(self.dialog)
        btype = self.filter_btype.currentText()
        #Group Box1 for General Parameeters
        group_box1 = QGroupBox()
        gb1 = QFormLayout()
        
        SamplingRate = QLineEdit(str(SR))
        
        gb1.addRow("Sampling Rate",SamplingRate)
        group_box1.setLayout(gb1)
        
        
        #Run Button
        Run = QPushButton('Run')
        Run.clicked.connect(plot_1)
        
        vbox.addWidget(group_box1)
        vbox.addWidget(Run)
        
        self.setLayout(vbox)
        self.dialog.setWindowTitle("Preprocess Settings")
        self.dialog.setWindowModality(Qt.ApplicationModal)
        self.dialog.exec_()
        
        
        
        
    def update_2(self):
        def plot_2():
            global SR, plot2_data, if_log
            SR = int(SamplingRate.text())
            LowFrequencyLimit = float(Low_freq.text())
            HighFrequencyLimit = float(High_freq.text())
            FAR = float(FrequencyAxisResolution.text())
            hop = int(tDS.text())
            WL = int(WindowLength.text())
            nW_InConceFT = int(NoWindowsInConceFT.text()) 
            WB = int(WindowBandwidth.text()) 
            n_ConceFT = int(NoConceFT.text())
            ifSecond = second_order.isChecked()
            ifSmooth = Smooth.isChecked()
            ifHemi = Hemi.isChecked()
            if_log[2] = self.zlog_cb2.currentIndex()
            
            
            if self.en_cb_2.currentIndex() == 0:
                z,y,_,_,_ = ConceFT_sqSTFT_C(signal_todo, LowFrequencyLimit, HighFrequencyLimit, FAR, hop, WL, nW_InConceFT, WB, 0, 0, 0,0)
            elif self.en_cb_2.currentIndex() == 1:    
                if ConceFT_cb.isChecked():
                    _, _, _, z, y = ConceFT_sqSTFT_C(signal_todo, LowFrequencyLimit, HighFrequencyLimit, FAR, hop, WL, nW_InConceFT, WB, n_ConceFT, ifSecond, ifSmooth, ifHemi) 
                else:
                    _, _, z, _, y = ConceFT_sqSTFT_C(signal_todo, LowFrequencyLimit, HighFrequencyLimit, FAR, hop, WL, nW_InConceFT, WB, 0, ifSecond, 0, 0)
            elif self.en_cb_2.currentIndex() == 2:
                if ConceFT_cb.isChecked():
                    _, _, _, z, y = ConceFT_rsSTFT_C(signal_todo, LowFrequencyLimit, HighFrequencyLimit, FAR, hop, WL, nW_InConceFT, WB, n_ConceFT) 
                else:
                    _, _, z, _, y = ConceFT_rsSTFT_C(signal_todo, LowFrequencyLimit, HighFrequencyLimit, FAR, hop, WL, nW_InConceFT, WB, 0)     
            
            
            plot2_data = []
            plot2_data.append(t_todo),plot2_data.append(y),plot2_data.append(z)
    
            imageSQ(self.figure, self.ax2, self.colorbar_2, t_todo, y*SR, np.abs(z), 99.5,if_log[1], 'Greys')
            self.ax2_signal.clear()
            self.ax2_signal.plot(t_todo,signal_todo, '-', color = 'r', linewidth = 0.3)
            self.ax2.set_title('STFT')
            self.ax2.set_xlabel('Time')
            self.ax2.set_ylabel('Frequency')
            self.ax2.grid(linestyle='dotted')
            self.sl_valuechange()
            self.dialog.accept()
            
            
            
        def ConceFT_cb_state(b):
                if b.isChecked() == True:
                    group_box4.setEnabled(True)
                else:
                    group_box4.setEnabled(False)
        
        
        
        
        #Confirm Data Loaded
        if signal_todo == []:
            print ('data not loaded')
        
            
        #Create Dialog for settings
        self.dialog = QDialog()
        vbox = QVBoxLayout(self.dialog)
        
        #Group Box1 for General Parameters
        group_box1 = QGroupBox()
        gb1 = QFormLayout()
        SamplingRate = QLineEdit(str(SR))
        SamplingRate.setValidator(QIntValidator())
        WindowLength = QLineEdit('377')
        WindowLength.setValidator(QIntValidator())
        WindowBandwidth = QLineEdit('10') 
        WindowBandwidth.setValidator(QIntValidator())
        tDS = QLineEdit('1')
        tDS.setValidator(QIntValidator())
        
        gb1.addRow("Sampling Rate",SamplingRate)
        gb1.addRow("Window Length",WindowLength)
        gb1.addRow("Window Bandwidth",WindowBandwidth)
        gb1.addRow("tDS",tDS)
        
        group_box1.setLayout(gb1)
        
       
       
        #SQ/RS box
        group_box2 = QGroupBox()
        gb2 = QFormLayout()
        FrequencyAxisResolution = QLineEdit('0.001')
        FrequencyAxisResolution.setValidator(QDoubleValidator())
        Low_freq = QLineEdit('0')
        Low_freq.setValidator(QDoubleValidator())
        High_freq = QLineEdit('0.5')
        High_freq.setValidator(QDoubleValidator())
        Smooth = QCheckBox()
        Hemi = QCheckBox()
        gb2.addRow("Frequency Axis Resolution", FrequencyAxisResolution)
        gb2.addRow("Low Frequency Limit",Low_freq)
        gb2.addRow("High_Frequency Limit",High_freq)
        gb2.addRow("Smooth",Smooth)
        gb2.addRow("Hemi",Hemi)
        group_box2.setLayout(gb2)
        
        
        #Order Box
        group_box3 = QGroupBox("&Order")
        first_order = QRadioButton("1st")
        second_order = QRadioButton("2nd")
        first_order.setChecked(True)
        gb3 = QVBoxLayout()
        gb3.addWidget(first_order)
        gb3.addWidget(second_order)
        group_box3.setLayout(gb3)
        
        
        #ConceFT Confirm
        ConceFT_cb = QCheckBox("Multitaper")
        ConceFT_cb.stateChanged.connect(lambda:ConceFT_cb_state(ConceFT_cb))
        group_box4 = QGroupBox()
        gb4 = QFormLayout()
        NoWindowsInConceFT = QLineEdit('2')
        NoWindowsInConceFT.setValidator(QIntValidator())
        NoConceFT = QLineEdit('20')
        NoConceFT.setValidator(QIntValidator())
        gb4.addRow("Number of windows in ConceFT", NoWindowsInConceFT)
        gb4.addRow("Number of ConceFT", NoConceFT)
        group_box4.setLayout(gb4)
        group_box4.setEnabled(False)
        
        #Run Button
        Run = QPushButton('Run')
        Run.clicked.connect(plot_2)
        
        #Add Widgets and Set Layout
        vbox.addWidget(group_box1)
        if self.en_cb_2.currentIndex() != 0:
            vbox.addWidget(group_box2)
            vbox.addWidget(ConceFT_cb)
            vbox.addWidget(group_box4)
        if self.en_cb_2.currentIndex() == 1:
            vbox.addWidget(group_box3)
        vbox.addWidget(Run)
        
        #Exececute
        
        self.setLayout(vbox)
        self.dialog.setWindowTitle("STFT Settings")
        self.dialog.setWindowModality(Qt.ApplicationModal)
        self.dialog.exec_()
        
        
    def update_3(self):
        def plot_3():
            global SR, plot3_data, if_log
            #Values for parameters
            SR = int(SamplingRate.text())
            LowFrequencyLimit = float(Low_freq.text())
            HighFrequencyLimit = float(High_freq.text())
            FAR = float(FrequencyAxisResolution.text())
            n_ConceFT = int(NoConceFT.text())
            ifSmooth = Smooth.isChecked()
            ifHemi = Hemi.isChecked()
            if_log[3] = self.zlog_cb3.currentIndex()
            #Valuese for opts
            opts_in.motherwavelet = MW_type
            opts_in.CENTER = float(CENTER.text())
            opts_in.FWHM = float(FWHM.text())
            opts_in.beta = float(beta.text())
            opts_in.gam = float(gam.text())
            opts_in.dim = int(dim.text())
            if MW_type == 'morse':
                opts_in.k = int(k.text())
            elif MW_type == 'morse-a':
                opts_in.k = k_morse_a.currentIndex()
            
            
            if self.en_cb_3.currentIndex() == 0:
                z, _, _, _, y = ConceFT_CWT(t_todo, signal_todo, LowFrequencyLimit, HighFrequencyLimit, FAR, 0, opts_in, 0, 0) 
            elif self.en_cb_3.currentIndex() == 1:    
                if ConceFT_cb.isChecked():
                    _, _, _, z, y = ConceFT_CWT(t_todo, signal_todo, LowFrequencyLimit, HighFrequencyLimit, FAR, n_ConceFT, opts_in, ifSmooth, ifHemi)  
                else:
                    _, _, z, _, y = ConceFT_CWT(t_todo, signal_todo, LowFrequencyLimit, HighFrequencyLimit, FAR, 0, opts_in, 0, 0) 
            
            plot3_data = []
            plot3_data.append(t_todo),plot3_data.append(y),plot3_data.append(z)
            self.colorbar_3 = imageSQ(self.figure, self.ax3, self.colorbar_3, t_todo, y*SR, np.abs(z), 99.5,if_log[2], 'Greys')
            self.ax3.set_title('CWT')
            self.ax3.set_xlabel('Time')
            self.ax3.set_ylabel('Frequency')
            self.ax3.grid(linestyle='dotted')
            self.sl_valuechange()
            self.dialog.accept()
            self.canvas.draw()
            
            
        def ConceFT_cb_state(b):
                if b.isChecked() == True:
                    group_box4.setEnabled(True)
                else:
                    group_box4.setEnabled(False)
        
        
        class opts_in:
            motherwavelet = 'Cinfc'
            CENTER = 1
            FWHM = 0.3
        
        if signal_todo == []:
            print ('data not loaded')
            return
        
        #Create Dialog for Settings
        self.dialog = QDialog()
        vbox = QVBoxLayout(self.dialog)
        
        #Group Box1 for General Parameeters
        group_box1 = QGroupBox()
        gb1 = QFormLayout()
        SamplingRate = QLineEdit(str(SR))
        SamplingRate.setValidator(QIntValidator())
        gb1.addRow("Sampling Rate",SamplingRate)
        group_box1.setLayout(gb1)
        
        
        #Multitaper box
        group_box2 = QGroupBox()
        gb2 = QFormLayout()
        FrequencyAxisResolution = QLineEdit('0.001')
        FrequencyAxisResolution.setValidator(QDoubleValidator())
        Low_freq = QLineEdit('0')
        Low_freq.setValidator(QDoubleValidator())
        High_freq = QLineEdit('0.5')
        High_freq.setValidator(QDoubleValidator())
        Smooth = QCheckBox()
        Hemi = QCheckBox()
        gb2.addRow("Frequency Axis Resolution", FrequencyAxisResolution)
        gb2.addRow("Low Frequency Limit",Low_freq)
        gb2.addRow("High_Frequency Limit",High_freq)
        gb2.addRow("Smooth",Smooth)
        gb2.addRow("Hemi",Hemi)
        group_box2.setLayout(gb2)
        
        
        #ConceFT Confirm
        ConceFT_cb = QCheckBox("Multitaper")
        ConceFT_cb.stateChanged.connect(lambda:ConceFT_cb_state(ConceFT_cb))
        group_box4 = QGroupBox()
        gb4 = QFormLayout()
        NoConceFT = QLineEdit('20')
        NoConceFT.setValidator(QIntValidator())
        gb4.addRow("Number of ConceFT", NoConceFT)
        group_box4.setLayout(gb4)
        group_box4.setEnabled(False)
        
        
        #Mother Wavelet Group Box
        MW_type = self.MW_cb3.currentText()
        MW_group_box = QGroupBox('Mother Wavelet')
        
        MW_gb = QFormLayout()
        l1 = QLabel()
        l1.setText(MW_type)
        CENTER = QLineEdit('1')
        CENTER.setValidator(QDoubleValidator())
        FWHM = QLineEdit('0.3')
        FWHM.setValidator(QDoubleValidator())
        dim = QLineEdit('1')
        dim.setValidator(QIntValidator())
        beta = QLineEdit('1')
        beta.setValidator(QDoubleValidator())
        gam = QLineEdit('1')
        gam.setValidator(QDoubleValidator())
        k = QLineEdit('0')
        k.setValidator(QIntValidator())
        k_morse_a = QComboBox()
        k_morse_a.addItems(['0','1'])
        
        MW_gb.addWidget(l1)
        if (MW_type == 'Cinfc') or (MW_type == 'gaussian'):
            MW_gb.addRow('Center', CENTER)
            MW_gb.addRow('FWHM',FWHM)
        elif MW_type in ['morse', 'morse-a', 'morse-b','morse-c'] :
            MW_gb.addRow('beta',beta)
            MW_gb.addRow('gam',gam)
            if MW_type in ['morse-b', 'morse-c']:
                MW_gb.addRow('dim',dim)
            elif MW_type == 'morse':
                MW_gb.addRow('k',k)
            elif MW_type == 'morse-a':
                MW_gb.addRow('k',k_morse_a)
        MW_group_box.setLayout(MW_gb)
        
        #Run Button
        Run = QPushButton('Run')
        Run.clicked.connect(plot_3)
        
        vbox.addWidget(MW_group_box)
        vbox.addWidget(group_box1)
        if self.en_cb_3.currentIndex() == 1:
            vbox.addWidget(group_box2)
            vbox.addWidget(ConceFT_cb)
            vbox.addWidget(group_box4)
        vbox.addWidget(Run)
        
        
        self.setLayout(vbox)
        self.dialog.setWindowTitle("CWT Settings")
        self.dialog.setWindowModality(Qt.ApplicationModal)
        self.dialog.exec_()
        
    def update_4(self):
        def plot_4():
            global SR, plot4_data, if_log
            #Values for parameters
            SR = int(SamplingRate.text())
            LowFrequencyLimit = float(Low_freq.text())
            HighFrequencyLimit = float(High_freq.text())
            FAR = float(FrequencyAxisResolution.text())
            ifSmooth = Smooth.isChecked()
            ifHemi = Hemi.isChecked()
            if_log[4] = self.zlog_cb4.currentIndex()
            if self.en_cb_4.currentIndex()==0:
                z,y,_,_ = stran(signal_todo,t_todo,0)
            if self.en_cb_4.currentIndex()==1:
                _,_,z,y = stran(signal_todo,t_todo,1,LowFrequencyLimit,HighFrequencyLimit,FAR,ifSmooth,ifHemi)
            
            plot4_data = []
            plot4_data.append(t_todo),plot4_data.append(y),plot4_data.append(z)
            self.colorbar_4 = imageSQ(self.figure, self.ax4, self.colorbar_4, t_todo, y*SR, np.abs(z), 99.5,if_log[3], 'Greys')
            self.ax4.set_title('S Transform')
            self.ax4.set_xlabel('Time')
            self.ax4.set_ylabel('Frequency')
            self.ax4.grid(linestyle='dotted')
            self.sl_valuechange()
            self.dialog.accept()
            self.canvas.draw()
            
            
        
        if signal_todo == []:
            print ('data not loaded')
            return
        
        #Create Dialog for Settings
        self.dialog = QDialog()
        vbox = QVBoxLayout(self.dialog)
        
        #Group Box1 for General Parameeters
        group_box1 = QGroupBox()
        gb1 = QFormLayout()
        SamplingRate = QLineEdit(str(SR))
        SamplingRate.setValidator(QIntValidator())
        gb1.addRow("Sampling Rate",SamplingRate)
        group_box1.setLayout(gb1)
        
        
        #SQ box
        group_box2 = QGroupBox()
        gb2 = QFormLayout()
        FrequencyAxisResolution = QLineEdit('0.001')
        FrequencyAxisResolution.setValidator(QDoubleValidator())
        Low_freq = QLineEdit('0')
        Low_freq.setValidator(QDoubleValidator())
        High_freq = QLineEdit('0.5')
        High_freq.setValidator(QDoubleValidator())
        Smooth = QCheckBox()
        Hemi = QCheckBox()
        gb2.addRow("Frequency Axis Resolution", FrequencyAxisResolution)
        gb2.addRow("Low Frequency Limit",Low_freq)
        gb2.addRow("High_Frequency Limit",High_freq)
        gb2.addRow("Smooth",Smooth)
        gb2.addRow("Hemi",Hemi)
        group_box2.setLayout(gb2)
        #Run Button
        Run = QPushButton('Run')
        Run.clicked.connect(plot_4)
        
        vbox.addWidget(group_box1)
        if self.en_cb_4.currentIndex() == 1:
            vbox.addWidget(group_box2)
        vbox.addWidget(Run)
        
        
        self.setLayout(vbox)
        self.dialog.setWindowTitle("S Transform Settings")
        self.dialog.setWindowModality(Qt.ApplicationModal)
        self.dialog.exec_()


           
    def update_5(self):
        def plot_5():
            global SR, plot5_data,if_log
            #Values for parameters
            SR = int(SamplingRate.text())
            FAR = float(FrequencyAxisResolution.text())
            LowFrequencyLimit = float(Low_freq.text())
            HighFrequencyLimit = float(High_freq.text())
            ifSmooth = Smooth.isChecked()
            ifHemi = Hemi.isChecked()
            if_log[5]=self.zlog_cb5.currentIndex()
            if self.en_cb_5.currentIndex()==0:
                z,y,_,_ = wvd(signal_todo,None,0)
            if self.en_cb_5.currentIndex()==1:
                _,_,z,y = wvd(signal_todo,None,1,FAR, 0, True, LowFrequencyLimit,HighFrequencyLimit,ifSmooth,ifHemi)
            
            plot5_data = []
            plot5_data.append(t_todo),plot5_data.append(y),plot5_data.append(z)
            self.colorbar_5 = imageSQ(self.figure, self.ax5, self.colorbar_5, t_todo, y*SR, np.abs(z), 99.5,if_log[4], 'Greys')
            self.ax5.set_title('WVD')
            self.ax5.set_xlabel('Time')
            self.ax5.set_ylabel('Frequency')
            self.ax5.grid(linestyle='dotted')
            self.sl_valuechange()
            self.dialog.accept()
            self.canvas.draw()
        
        
        if signal_todo == []:
            print ('data not loaded')
            return
        
        #Create Dialog for Settings
        self.dialog = QDialog()
        vbox = QVBoxLayout(self.dialog)
        
        #Group Box1 for General Parameeters
        group_box1 = QGroupBox()
        gb1 = QFormLayout()
        SamplingRate = QLineEdit(str(SR))
        SamplingRate.setValidator(QIntValidator())
        gb1.addRow("Sampling Rate",SamplingRate)
        group_box1.setLayout(gb1)
        
        #sq box
        group_box2 = QGroupBox()
        gb2 = QFormLayout()
        FrequencyAxisResolution = QLineEdit('0.001')
        FrequencyAxisResolution.setValidator(QDoubleValidator())
        Low_freq = QLineEdit('0')
        Low_freq.setValidator(QDoubleValidator())
        High_freq = QLineEdit('0.5')
        High_freq.setValidator(QDoubleValidator())
        Smooth = QCheckBox()
        Hemi = QCheckBox()
        gb2.addRow("Frequency Axis Resolution", FrequencyAxisResolution)
        gb2.addRow("Low Frequency Limit",Low_freq)
        gb2.addRow("High_Frequency Limit",High_freq)
        gb2.addRow("Smooth",Smooth)
        gb2.addRow("Hemi",Hemi)
        group_box2.setLayout(gb2)
        
        #Run Button
        Run = QPushButton('Run')
        Run.clicked.connect(plot_5)
        
        vbox.addWidget(group_box1)
        
        if self.en_cb_5.currentIndex() == 1:
            vbox.addWidget(group_box2)
        
        vbox.addWidget(Run)
        
        
        self.setLayout(vbox)
        self.dialog.setWindowTitle("WVD Settings")
        self.dialog.setWindowModality(Qt.ApplicationModal)
        self.dialog.exec_()
    
        
        



def main():
   app = QApplication(sys.argv)
   ex = Window()
   ex.showMaximized()
   sys.exit(app.exec_())
	
if __name__ == '__main__':
   main()