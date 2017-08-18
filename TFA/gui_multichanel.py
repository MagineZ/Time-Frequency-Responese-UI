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
import os
from os.path import basename
from copy import deepcopy
from scipy import signal

signal_dic = {};t_min_dic = {};t_max_dic={}
plot_count = 0
total_plot = 10
plot_row_n = 1
view_wl = 20
t_min = 0
t_max = 99
max_order = 0
fig_len = 19
fig_wid = 9


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
            global signal_dic,t_min_dic,t_max_dic
            self.name = QFileDialog.getOpenFileName(self, 'Open Data','' ,"*.txt *.csv")
            
            if self.name[0]:
                data = np.genfromtxt(self.name[0], delimiter=',')
                n = len(data.shape)
                if n == 1:
                    if basename(self.name[0]) not in signal_dic:
                        self.signal_list.addItem(basename(self.name[0]))
                    signal_dic[basename(self.name[0])] = [np.arange(0,len(data)),np.arange(0,len(data)),data, data,[False]*10,[False]*10,[0]*10,[[]]*10]
                    t_min_dic[basename(self.name[0])] = 0; t_max_dic[basename(self.name[0])] = len(data)-1
                    update_time_lim()
                    
                elif n==2:
                    if basename(self.name[0]) not in signal_dic:
                        self.signal_list.addItem(basename(self.name[0]))
                    signal_dic[basename(self.name[0])] = [data[:,0], data[:,0], data[:,1],data[:,1],[False]*10,[False]*10,[0]*10,[[]]*10]
                    t_min_dic[basename(self.name[0])] = np.min( data[:,0]) ; t_max_dic[basename(self.name[0])] = np.max( data[:,0])
                    update_time_lim()
                else:
                    print ('Wrong Format')
        
        elif q.text() == "Plot Selection":
            global signal_list
            self.ps_dialog = QDialog()
            self.plot_cb = []
            dialog_box = QVBoxLayout(self.ps_dialog)
            for i in range(0,total_plot):
                self.plot_cb.append(QCheckBox())
            self.signal_list = QListWidget()
            self.signal_list.itemClicked.connect(self.signal_list_itemClicked)
            for i in signal_dic:
                self.signal_list.addItem(i)
            self.plot_cb[0].setText('Signal View')
            self.plot_cb[1].setText('Power Spectrum')
            self.plot_cb[2].setText('STFT')
            self.plot_cb[3].setText('CWT')
            self.plot_cb[4].setText('S Transform')
            self.plot_cb[5].setText('WVD')
            self.plot_cb[6].setText('Ceps')
            self.plot_cb[7].setText('TCeps')
            self.plot_cb[8].setText('TFRR')
            self.plot_cb[9].setText('RTFR')
            
            
            
            vbox = QGridLayout()
            self.plot_cb_title = QLabel('Plot Selection')
            vbox.addWidget(self.plot_cb_title,1,1)
            self.plot_cb[0].stateChanged.connect(lambda:self.btnstate(self.plot_cb[0]))
            self.plot_cb[1].stateChanged.connect(lambda:self.btnstate(self.plot_cb[1]))
            self.plot_cb[2].stateChanged.connect(lambda:self.btnstate(self.plot_cb[2]))
            self.plot_cb[3].stateChanged.connect(lambda:self.btnstate(self.plot_cb[3]))
            self.plot_cb[4].stateChanged.connect(lambda:self.btnstate(self.plot_cb[4]))
            self.plot_cb[5].stateChanged.connect(lambda:self.btnstate(self.plot_cb[5]))
            self.plot_cb[6].stateChanged.connect(lambda:self.btnstate(self.plot_cb[6]))
            self.plot_cb[7].stateChanged.connect(lambda:self.btnstate(self.plot_cb[7]))
            self.plot_cb[8].stateChanged.connect(lambda:self.btnstate(self.plot_cb[8]))
            self.plot_cb[9].stateChanged.connect(lambda:self.btnstate(self.plot_cb[9]))
            
            for i in range(0,total_plot):
                vbox.addWidget(self.plot_cb[i],i+2,1)
                
            
    
                
            
            self.plot_order_title = QLabel('Plot Order')
            self.plot_order = []
            vbox.addWidget(QLabel('Plot Order'),1,2)
            for i in range(0,total_plot):
                self.plot_order.append(QLineEdit('0'))
                self.plot_order[i].setValidator(QIntValidator())
                vbox.addWidget(self.plot_order[i],i+2,2)
                
            self.plot_order[0].textChanged.connect(lambda:self.plot_order_changed(self.plot_order[0]))
            self.plot_order[1].textChanged.connect(lambda:self.plot_order_changed(self.plot_order[1]))
            self.plot_order[2].textChanged.connect(lambda:self.plot_order_changed(self.plot_order[2]))
            self.plot_order[3].textChanged.connect(lambda:self.plot_order_changed(self.plot_order[3]))
            self.plot_order[4].textChanged.connect(lambda:self.plot_order_changed(self.plot_order[4]))
            self.plot_order[5].textChanged.connect(lambda:self.plot_order_changed(self.plot_order[5]))
            self.plot_order[6].textChanged.connect(lambda:self.plot_order_changed(self.plot_order[6]))
            self.plot_order[7].textChanged.connect(lambda:self.plot_order_changed(self.plot_order[7]))
            self.plot_order[8].textChanged.connect(lambda:self.plot_order_changed(self.plot_order[8]))
            self.plot_order[9].textChanged.connect(lambda:self.plot_order_changed(self.plot_order[9]))
            
            self.PlotRow = QLineEdit(str(plot_row_n))
            self.PlotRow.setValidator(QIntValidator())
            self.PlotRow_title = QLabel('Number of Rows in Plot')
            vbox.addWidget(self.PlotRow_title,total_plot+2,1)
            vbox.addWidget(self.PlotRow,total_plot+2,2)
            
            self.fig_len = QLineEdit(str(fig_len)); self.fig_len.setValidator(QIntValidator())
            self.fig_wid = QLineEdit(str(fig_wid)); self.fig_wid.setValidator(QIntValidator())
            self.fig_size_box = QHBoxLayout()
            self.fig_size_box.addWidget(self.fig_wid);self.fig_size_box.addWidget(self.fig_len)
            self.fig_size_box_title = QLabel('Figure Size')
            vbox.addWidget(self.fig_size_box_title,total_plot+3,1)
            vbox.addLayout(self.fig_size_box,total_plot+3,2)
            
            hbox = QHBoxLayout()
            Add_signal = QPushButton('Add')
            Delete_signal = QPushButton('Delete')
            Done = QPushButton('Done')
            Add_signal.clicked.connect(self.load_data)
            Delete_signal.clicked.connect(self.delete_data)
            Done.clicked.connect(self.plot_selection_Done)
            hbox.addWidget(Add_signal);hbox.addWidget(Delete_signal);hbox.addWidget(Done)
            
            dialog_box.addWidget(self.signal_list)
            dialog_box.addLayout(vbox)
            dialog_box.addLayout(hbox)
        
            self.setLayout(dialog_box)
            self.ps_dialog.setWindowTitle("Plot Selection")
            self.ps_dialog.setWindowModality(Qt.ApplicationModal)
            self.ps_dialog.exec_()
        
    
    def signal_list_itemClicked(self):
        global signal_dic
        items = self.signal_list.selectedItems()
        for i in items:
            for j in range(0,total_plot):
                self.plot_cb[j].setChecked(signal_dic[i.text()][4][j])
                self.plot_order[j].setText(str(signal_dic[i.text()][6][j]))
        
        
    def plot_order_changed(self,b):
        global signal_dic
        items = self.signal_list.selectedItems()
        for i in items:
            if b.text() != '':
                signal_dic[i.text()][6][self.plot_order.index(b)]=int(b.text())
            print (signal_dic[i.text()][6])
    
    
    def btnstate(self,b):
        global signal_dic
        items = self.signal_list.selectedItems()
        if b.isChecked():
            for i in items: 
                signal_dic[i.text()][4][self.plot_cb.index(b)]=True
                print(i.text()+str(signal_dic[i.text()][4]))
        else:
            for i in items: 
                signal_dic[i.text()][4][self.plot_cb.index(b)]=False 
                print(i.text()+str(signal_dic[i.text()][4]))
        return
        
     
    
    def load_data(self):
        global signal_dic, t_min, t_max
        self.name = QFileDialog.getOpenFileName(self, 'Open Data','' ,"*.txt *.csv")
        if self.name[0]:
            data = np.genfromtxt(self.name[0], delimiter=',')
            n = len(data.shape)
            if n == 1:
                if basename(self.name[0]) not in signal_dic:
                    self.signal_list.addItem(basename(self.name[0]))
                signal_dic[basename(self.name[0])] = [np.arange(0,len(data)),np.arange(0,len(data)),data, data,[False]*10,[False]*10,[0]*10,[[]]*10]
                t_min_dic[basename(self.name[0])] = 0; t_max_dic[basename(self.name[0])] = len(data)-1
                update_time_lim()
            elif n==2:
                if basename(self.name[0]) not in signal_dic:
                    self.signal_list.addItem(basename(self.name[0]))
                signal_dic[basename(self.name[0])] = [data[:,0], data[:,0], data[:,1],data[:,1],[False]*10,[False]*10,[0]*10,[[]]*10]
                t_min_dic[basename(self.name[0])] = np.min( data[:,0]) ; t_max_dic[basename(self.name[0])] = np.max( data[:,0])
                update_time_lim()
            else:
                print ('Wrong Format')
        
        
        
        
    def delete_data(self):
        global signal_dic,t_min_dic,t_max_dic
        items = self.signal_list.selectedItems()
        for i in items:
            del signal_dic[i.text()]
            del t_min_dic[i.text()]
            del t_max_dic[i.text()]
            update_time_lim()
            self.signal_list.takeItem(self.signal_list.row(i))
        
        
        
    def plot_selection_Done(self):
        global plot_count, max_order, plot_row_n, fig_len, fig_wid
        plot_count = 0
        max_order = 0
        fig_len = int(self.fig_len.text())
        fig_wid = int(self.fig_wid.text())
        for i in signal_dic:
            plot_count+=sum(signal_dic[i][4])
            temp_max_order = max(signal_dic[i][6])
            if temp_max_order>max_order:
                max_order = temp_max_order 
        plot_row_n = int(self.PlotRow.text())
        self.form_widget.run_plot()
        self.ps_dialog.accept()






# Central Widget
class FormWidget(QWidget):

    def __init__(self, parent):        
        super(FormWidget, self).__init__(parent)
        #Figure1
        
        self.figure_hide = plt.figure()
        self.mw_layout_vbox = QVBoxLayout()
        
        
        
        self.e_t = np.arange(0,100)
        self.e_y = np.arange(0,100)
        self.e_M = np.zeros((100,100))
        self.e_M[0,0] = 10           
        #self.e_cax = self.ax_hide.pcolorfast(self.e_t, self.e_y, self.e_M, cmap = 'Greys')
        
        
        #slider
        self.sl = QSlider(Qt.Horizontal)
        self.sl_min = QLabel('0')
        self.sl.setMinimum(0)
        self.sl_max = QLabel(str(100-view_wl))
        self.sl.setMaximum(100-view_wl)
        self.sl.setValue(0)
        self.sl_cur_val = QLabel(str(self.sl.value()))
        self.sl.setTickPosition(QSlider.TicksBelow)
        self.sl.valueChanged.connect(self.sl_valuechange)
        self.sl_wl = QLineEdit(str(view_wl))
        self.sl_wl.setFixedWidth(50)
        self.sl_wl.setValidator(QIntValidator())
        self.sl_wl.returnPressed.connect(self.sl_valuechange)
        
    
  
    def run_plot(self):
        self.clearLayout(self.mw_layout_vbox)
        self.signal_list = QListWidget()
        for i in signal_dic:
            self.signal_list.addItem(i)
        self.mw_layout_hbox = QHBoxLayout()
        self.figure = plt.figure()
        self.figure.set_size_inches(fig_len,fig_wid)
        self.canvas = FigureCanvas(self.figure)
        self.scrollarea = QScrollArea(self)
        #if max_order<7:
        #    self.scrollarea.setWidgetResizable(True)
        self.scrollarea.setWidget(self.canvas)
        
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.panel_base = QWidget()
        self.panel = QVBoxLayout()
        if max_order>3 and plot_row_n<2:
            self.nrow = 2
        else:
            self.nrow = plot_row_n
        self.ncol = int(np.ceil(max_order/self.nrow))
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
        self.lines = {}; self.lines_super = {}; self.SR = {}
        self.axes_hide = {}
        hide_count = 1
        for i in signal_dic:
            self.lines[i] = [[]]*total_plot
            self.lines_super[i]=[[]]*total_plot
            self.SR[i] = 1
            self.axes_hide[i] = self.figure_hide.add_subplot(len(signal_dic),1,hide_count)
            self.axes_hide[i].set_xlim([0,0+view_wl])
            hide_count+=1
        self.axes_record = [[]]*max_order
        self.axes_super = [[]]*max_order
        self.colorbars = [[]]*max_order
        self.if_log = [False]*max_order
        self.if_superimpose = [False]*max_order
        self.add_panels()
        
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
        #self.sl_panel_base.setFixedHeight(50)
        self.sl_panel_base.setLayout(self.sl_vlayout)
        self.panel_base.setLayout(self.panel)
        self.panel_base.setFixedWidth(150)

        self.scroll_panel = QScrollArea(self)
        self.scroll_panel.setWidget(self.panel_base)
        self.scroll_panel.setFixedWidth(150)
        # set the layout
        #self.figure.tight_layout()
        "self.figure.subplotpars.update(left=0.05, bottom=0.05, right=0.95, top=0.95, wspace=0.2, hspace=0.2)"
        
        
        self.mw_layout_vbox.addWidget(self.toolbar)
        self.mw_layout_hbox.addWidget(self.scroll_panel)
        self.figure_vbox = QVBoxLayout()
        self.figure_vbox.addWidget(self.scrollarea)
        self.figure_vbox.addWidget(self.sl_panel_base)
        self.mw_layout_hbox.addLayout(self.figure_vbox)
        self.mw_layout_vbox.addLayout(self.mw_layout_hbox)
        self.setLayout(self.mw_layout_vbox)
        self.sl_valuechange()
        
    
    def add_panels(self):
        if_add_panel = [False]*total_plot
        for i in signal_dic:
            for j in range(0,total_plot):
                if_add_panel[j] = if_add_panel[j] | signal_dic[i][4][j]
            if signal_dic[i][4][0] and signal_dic[i][6][0]>0:
                order = signal_dic[i][6][0]
                if not self.axes_record[order-1]:
                    self.axes_record[order-1] = self.figure.add_subplot(self.nrow,self.ncol,order,sharex = self.axes_hide[i])
                    self.axes_record[order-1].grid(linestyle='dotted')
                
                #if signal_dic[i][7][0]:
                self.axes_record[order-1].plot(signal_dic[i][1],signal_dic[i][3], '-', linewidth = 1,label=i)
                self.lines[i][0] = self.axes_record[order-1].lines[-1]
                #self.axes[i][0].set_xlim([np.min(plot1_data[0]),np.max(plot1_data[0])])
                self.axes_record[order-1].legend(shadow=True, fancybox=True)
                self.axes_record[order-1].set_title('Signal View')
                self.axes_record[order-1].set_xlabel('Time')
                if not self.colorbars[order-1]:
                    self.colorbars[order-1] = self.figure.colorbar(self.axes_hide[i].pcolorfast(self.e_t, self.e_y, self.e_M, cmap = 'Greys'), ax=self.axes_record[order-1], extend='max')
                    self.colorbars[order-1].remove()
            
            if signal_dic[i][4][1] and signal_dic[i][6][1]>0:
                order = signal_dic[i][6][1]
                if not self.axes_record[order-1]:
                    self.axes_record[order-1] = self.figure.add_subplot(self.nrow,self.ncol,order)
                    self.axes_record[order-1].grid(linestyle='dotted')
                
                self.axes_record[order-1].legend(shadow=True, fancybox=True)
                self.axes_record[order-1].set_title('Power Spectrum')
                self.axes_record[order-1].set_xlabel('Frequency')
                self.axes_record[order-1].set_ylabel('V**2')
                if not self.colorbars[order-1]:
                    self.colorbars[order-1] = self.figure.colorbar(self.axes_hide[i].pcolorfast(self.e_t, self.e_y, self.e_M, cmap = 'Greys'), ax=self.axes_record[order-1], extend='max')
                    self.colorbars[order-1].remove()
                    
            for k in range(2,total_plot):
                if signal_dic[i][4][k] and signal_dic[i][6][k]>0:
                    order = signal_dic[i][6][k]
                    if not self.axes_record[order-1]:
                        self.axes_record[order-1] = self.figure.add_subplot(self.nrow,self.ncol,order,sharex = self.axes_hide[i],sharey=self.axes_hide[i])
                    if not self.colorbars[order-1]:
                        self.colorbars[order-1] = self.figure.colorbar(self.axes_record[order-1].pcolorfast(self.e_t, self.e_y, self.e_M, cmap = 'Greys'), ax=self.axes_record[order-1], extend='max')        
                    self.axes_record[order-1].grid(linestyle='dotted')  
                    self.axes_record[order-1].set_xlabel('Time')
                    if k!= 6:
                        self.axes_record[order-1].set_ylabel('Frequency')
                    else:
                        self.axes_record[order-1].set_ylabel('quefrency')
                    if k == 2:self.axes_record[order-1].set_title('STFT: '+i)
                    if k == 3:self.axes_record[order-1].set_title('CWT: '+i)
                    if k == 4:self.axes_record[order-1].set_title('S Transform: '+i)
                    if k == 5:self.axes_record[order-1].set_title('WVD: '+i)
                    if k == 6:self.axes_record[order-1].set_title('Ceps: '+i)
                    if k == 7:self.axes_record[order-1].set_title('tCeps: '+i)
                    if k == 8:self.axes_record[order-1].set_title('TFRR: '+i)
                    if k == 9:self.axes_record[order-1].set_title('RTFR: '+i)

            
                
                
                
                
                
                
                
                
                
        
        if if_add_panel[0]:
            self.gb_signal = QGroupBox('Signal View')
            self.para_layout1 = QVBoxLayout()
            self.button_1 = QPushButton('Update')
            self.button_1.clicked.connect(self.update_signal)
            self.para_layout1.addWidget(self.button_1)
            self.gb_signal.setLayout(self.para_layout1)
            self.panel.addWidget(self.gb_signal)
        
        if if_add_panel[1]:
            self.gb_ps = QGroupBox('Power Spectrum')
            self.para_layout0 = QVBoxLayout()
            self.button_0 = QPushButton('Update')
            self.button_0.clicked.connect(self.update_ps)
            self.para_layout0.addWidget(self.button_0)
            self.gb_ps.setLayout(self.para_layout0)
            self.panel.addWidget(self.gb_ps)
        
        if if_add_panel[2]:
            self.gb_STFT = QGroupBox('STFT')
            self.para_layout2 = QVBoxLayout()
            self.zlog_cb2 = QComboBox()
            self.zlog_cb2.addItems(["normal scale","log scale"])
            self.en_cb_2 = QComboBox()
            self.en_cb_2.addItems(["Normal", "Synchrosquezzing", "Reassigment"])
            self.button_2 = QPushButton('Update')
            self.button_2.clicked.connect(self.update_STFT)
            self.para_layout2.addWidget(self.zlog_cb2)
            self.para_layout2.addWidget(self.en_cb_2)
            self.para_layout2.addWidget(self.button_2)
            self.gb_STFT.setLayout(self.para_layout2)
            self.panel.addWidget(self.gb_STFT)
            
        if if_add_panel[3]:
            self.gb_CWT = QGroupBox('CWT')
            self.para_layout3 = QVBoxLayout()
            self.zlog_cb3 = QComboBox()
            self.zlog_cb3.addItems(["normal scale","log scale"])
            self.en_cb_3 = QComboBox()
            self.en_cb_3.addItems(["Normal", "Synchrosquezzing"])
            self.MW_cb3 = QComboBox()
            self.MW_cb3.addItems(['Cinfc','morse', 'morse-a', 'morse-b', 'morse-c', 'morlet', 'gaussian', 'meyer', 'BL3'])
            self.button_3 = QPushButton('Update')
            self.button_3.clicked.connect(self.update_CWT)
            self.para_layout3.addWidget(self.zlog_cb3)
            self.para_layout3.addWidget(self.en_cb_3)
            self.para_layout3.addWidget(self.MW_cb3)
            self.para_layout3.addWidget(self.button_3)
            self.gb_CWT.setLayout(self.para_layout3)
            self.panel.addWidget(self.gb_CWT)
        
        if if_add_panel[4]:
            self.gb_stran = QGroupBox('S Transform')
            self.para_layout4 = QVBoxLayout()
            self.zlog_cb4 = QComboBox()
            self.zlog_cb4.addItems(["normal scale","log scale"])
            self.en_cb_4 = QComboBox()
            self.en_cb_4.addItems(["Normal", "Synchrosquezzing"])
            self.button_4 = QPushButton('Update')
            self.button_4.clicked.connect(self.update_stran)
            self.para_layout4.addWidget(self.zlog_cb4)
            self.para_layout4.addWidget(self.en_cb_4)
            self.para_layout4.addWidget(self.button_4)
            self.gb_stran.setLayout(self.para_layout4)
            self.panel.addWidget(self.gb_stran)
        
        if if_add_panel[5]:
            self.gb_wvd = QGroupBox('WVD')
            self.para_layout5 = QVBoxLayout()
            self.zlog_cb5 = QComboBox()
            self.zlog_cb5.addItems(["normal scale","log scale"])
            self.en_cb_5 = QComboBox()
            self.en_cb_5.addItems(["Normal", "Synchrosquezzing"])
            self.button_5 = QPushButton('Update')
            self.button_5.clicked.connect(self.update_WVD)
            self.para_layout5.addWidget(self.zlog_cb5)
            self.para_layout5.addWidget(self.en_cb_5)
            self.para_layout5.addWidget(self.button_5)
            self.gb_wvd.setLayout(self.para_layout5)
            self.panel.addWidget(self.gb_wvd)
            
        if if_add_panel[6]:
            self.gb_ceps = QGroupBox('Cepstrum')
            self.para_layout6 = QVBoxLayout()
            self.zlog_cb6 = QComboBox()
            self.zlog_cb6.addItems(["normal scale","log scale"])
            self.button_6 = QPushButton('Update')
            self.button_6.clicked.connect(lambda:self.update_CFPH('CEPS'))
            self.para_layout6.addWidget(self.zlog_cb6)
            self.para_layout6.addWidget(self.button_6)
            self.gb_ceps.setLayout(self.para_layout6)
            self.panel.addWidget(self.gb_ceps)
        
        if if_add_panel[7]:
            self.gb_tceps = QGroupBox('tCeps')
            self.para_layout7 = QVBoxLayout()
            self.zlog_cb7 = QComboBox()
            self.zlog_cb7.addItems(["normal scale","log scale"])
            self.button_7 = QPushButton('Update')
            self.button_7.clicked.connect(lambda:self.update_CFPH('TCEPS'))
            self.para_layout7.addWidget(self.zlog_cb7)
            self.para_layout7.addWidget(self.button_7)
            self.gb_tceps.setLayout(self.para_layout7)
            self.panel.addWidget(self.gb_tceps)
        
        if if_add_panel[8]:
            self.gb_tfrr = QGroupBox('TFRR')
            self.para_layout8 = QVBoxLayout()
            self.zlog_cb8 = QComboBox()
            self.zlog_cb8.addItems(["normal scale","log scale"])
            self.button_8 = QPushButton('Update')
            self.button_8.clicked.connect(lambda:self.update_CFPH('TFRR'))
            self.para_layout8.addWidget(self.zlog_cb8)
            self.para_layout8.addWidget(self.button_8)
            self.gb_tfrr.setLayout(self.para_layout8)
            self.panel.addWidget(self.gb_tfrr)
        
        if if_add_panel[9]:
            self.gb_rtfr = QGroupBox('RTFR')
            self.para_layout9 = QVBoxLayout()
            self.zlog_cb9 = QComboBox()
            self.zlog_cb9.addItems(["normal scale","log scale"])
            self.button_9 = QPushButton('Update')
            self.button_9.clicked.connect(lambda:self.update_CFPH('RTFR'))
            self.para_layout9.addWidget(self.zlog_cb9)
            self.para_layout9.addWidget(self.button_9)
            self.gb_rtfr.setLayout(self.para_layout9)
            self.panel.addWidget(self.gb_rtfr)
            
            
        
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

            
    def sl_valuechange(self):
        global view_wl
        view_wl = int(self.sl_wl.text())
        t_len = t_max-t_min+1
        if view_wl>=(t_len):
            view_wl = int(t_len)
            self.sl.setMaximum(t_max+1)
            self.sl.setMinimum(t_min)
            self.sl_max.setText(str(t_max))
            self.sl_min.setText(str(t_min))
            self.sl_cur_val.setText(str(t_min))
            for i in signal_dic:    
                self.axes_hide[i].set_xlim([t_min,t_max])
            self.sl.setValue(t_min)
        else:
            pos = self.sl.value()
            if (pos<t_min) or (pos>t_max-view_wl):
                pos = t_min
                
            self.sl.setValue(pos)    
            self.sl_cur_val.setText(str(pos))
            for i in self.axes_hide:
                self.axes_hide[i].set_xlim([pos,pos+view_wl])
            self.sl.setMaximum(t_max-view_wl)
            self.sl.setMinimum(t_min)
            self.sl_max.setText(str(t_max-view_wl))
            self.sl_min.setText(str(t_min))
    
    
        
    def update_pre(self):
        def plot_pre():
            global signal_dic, t_min_dic, t_max_dic
            DT = detrend.isChecked()
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
            items = self.signal_list.selectedItems()
            for i in items:
                self.SR[i.text()] = int(SamplingRate.text())
                if btype!='None':
                    signal_dic[i.text()][1], signal_dic[i.text()][3] = signal_filter(signal_dic[i.text()][0],signal_dic[i.text()][2],self.SR[i.text()],DT,DS, btype,ftype.currentText(),N,cutoff,rp,rs)
                else:
                    signal_dic[i.text()][1], signal_dic[i.text()][3] = signal_filter(signal_dic[i.text()][0],signal_dic[i.text()][2],self.SR[i.text()],DT,DS, btype)
                t_min_dic[i.text()] = np.min(signal_dic[i.text()][1]) ; t_max_dic[i.text()] = np.max( signal_dic[i.text()][1])
                update_time_lim()
        
        def close_window():
            self.dialog.accept()
            
        if signal_dic == {}:
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
        SamplingRate = QLineEdit(str(100))
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
        self.signal_list = QListWidget()
        for i in signal_dic:
            self.signal_list.addItem(i)
        gb_pre.addWidget(self.signal_list)
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
        Run = QPushButton('Update')
        Run.clicked.connect(plot_pre)
        Close = QPushButton('Close')
        Close.clicked.connect(close_window)
        vbox.addWidget(group_box_pre)
        vbox.addWidget(Run)
        vbox.addWidget(Close)
        
        self.setLayout(vbox)
        self.dialog.setWindowTitle("Preprocess Settings")
        self.dialog.setWindowModality(Qt.ApplicationModal)
        self.dialog.exec_()
        
    
    def update_signal(self):
        def plot_signal(delete_line):
            global signal_dic
            items = self.signal_list.selectedItems()
            for j in items:
                i = j.text()
                SR = self.SR[i]
                order = signal_dic[i][6][0]
                #if signal_dic[i][7][0]:
                if self.lines[i][0]:   
                    self.lines[i][0].remove()
                    self.lines[i][0] = []
                if not delete_line:
                    self.axes_record[order-1].plot(signal_dic[i][1],signal_dic[i][3], '-', linewidth = 1,label=i)
                    self.lines[i][0] = self.axes_record[order-1].lines[-1]
                    self.axes_record[order-1].legend(shadow=True, fancybox=True)
                #self.axes[i][0].set_xlim([np.min(plot1_data[0]),np.max(plot1_data[0])])
                    
            self.sl_valuechange()
            self.canvas.draw()
            
        def close_window():
            self.dialog.accept()
        
        if not signal_dic:
            print ('data not loaded')
            return
        self.dialog = QDialog()
        self.signal_list = QListWidget()
        for i in signal_dic:
            self.signal_list.addItem(i)
        vbox = QVBoxLayout(self.dialog)
        #Group Box1 for General Parameeters
        group_box1 = QGroupBox()
        gb1 = QFormLayout()
        group_box1.setLayout(gb1)

        
        #Run Button
        Run = QPushButton('Run')
        Run.clicked.connect(lambda: plot_signal(0))
        Delete = QPushButton('Delete')
        Delete.clicked.connect(lambda:plot_signal(1))
        Close = QPushButton('Close')
        Close.clicked.connect(close_window)
        vbox.addWidget(self.signal_list)
        vbox.addWidget(group_box1)
        vbox.addWidget(Run)
        vbox.addWidget(Delete)
        vbox.addWidget(Close)
        
        self.setLayout(vbox)
        self.dialog.setWindowTitle("Signal View Settings")
        self.dialog.setWindowModality(Qt.ApplicationModal)
        self.dialog.exec_()
    
    def update_ps(self):
        def plot_ps(delete_line):
            global signal_dic
            items = self.signal_list.selectedItems()
            for j in items:
                i = j.text()
                SR = self.SR[i]
                WL = int(WindowLength.text())
                scale = Scaling.currentText()
                order = signal_dic[i][6][1]
                #if signal_dic[i][7][0]:
                if self.lines[i][1]:   
                    self.lines[i][1].remove()
                    self.lines[i][1] = []
                if not delete_line:
                    x, y = signal.periodogram(signal_dic[i][3], SR,nfft = WL, scaling = scale)
                    self.axes_record[order-1].semilogy(x[1:], y[1:],linewidth = 1,label = i)
                    self.axes_record[order-1].legend(shadow=True, fancybox=True)
                    self.lines[i][1] = self.axes_record[order-1].lines[-1]
                #self.axes[i][0].set_xlim([np.min(plot1_data[0]),np.max(plot1_data[0])])
                    
            self.sl_valuechange()
            self.canvas.draw()
            
        def close_window():
            self.dialog.accept()
        
        if not signal_dic:
            print ('data not loaded')
            return
        self.dialog = QDialog()
        self.signal_list = QListWidget()
        for i in signal_dic:
            self.signal_list.addItem(i)
        vbox = QVBoxLayout(self.dialog)
        #Group Box1 for General Parameeters
        group_box1 = QGroupBox()
        gb1 = QFormLayout()
        WindowLength = QLineEdit('377')
        WindowLength.setValidator(QIntValidator())
        Scaling = QComboBox()
        Scaling.addItems(['spectrum','density'])
        gb1.addRow('Window Length', WindowLength)
        gb1.addRow('Scaling', Scaling)
        group_box1.setLayout(gb1)

        
        #Run Button
        Run = QPushButton('Run')
        Run.clicked.connect(lambda:plot_ps(0))
        Delete = QPushButton('Delete')
        Delete.clicked.connect(lambda:plot_ps(1))
        Close = QPushButton('Close')
        Close.clicked.connect(close_window)
        vbox.addWidget(self.signal_list)
        vbox.addWidget(group_box1)
        vbox.addWidget(Run)
        vbox.addWidget(Delete)
        vbox.addWidget(Close)
        
        self.setLayout(vbox)
        self.dialog.setWindowTitle("Power Spectrum Settings")
        self.dialog.setWindowModality(Qt.ApplicationModal)
        self.dialog.exec_()
        
    def update_STFT(self):
        def plot_2():
            global signal_dic
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
            items = self.signal_list.selectedItems()
            for j in items:
                i = j.text(); order = signal_dic[i][6][2]
                self.if_superimpose[order-1] = SuperImpose.isChecked()
                SR = self.SR[i]; self.if_log[order-1] = self.zlog_cb2.currentIndex()
                signal_todo = signal_dic[i][3]
                t_todo = signal_dic[i][1]
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
                
                self.colorbars[order-1] = imageSQ(self.figure, self.axes_record[order-1], self.colorbars[order-1], t_todo, y*SR, np.abs(z), 99.5,self.if_log[order-1], 'Greys')
                self.axes_record[order-1].legend(shadow=True, fancybox=True)
        
            
                if self.if_superimpose[order-1]:
                    if not self.axes_super[order-1]:
                        self.axes_super[order-1] = self.axes_record[order-1].twinx()
                        self.axes_super[order-1].plot(signal_dic[i][1],signal_dic[i][3], '-', linewidth = 1,label=i,color='r')
                        self.lines_super[i][2] = self.axes_super[order-1].lines[-1]
                        #self.axes_super[order-1].set_title('Superimpose: '+i,loc = 'right')
                    else:
                        if self.lines_super[i][2]:
                            self.lines_super[i][2].remove()
                            self.lines_super[i][2] = []
                        self.axes_super[order-1].plot(signal_dic[i][1],signal_dic[i][3], '-', linewidth = 1,label=i,color='r')
                        self.lines_super[i][2] = self.axes_super[order-1].lines[-1]
                        #self.axes_super[order-1].set_title('Superimpose: '+i, loc = 'right')
                else:
                    if self.axes_super[order-1] and self.lines_super[i][2]:
                        self.lines_super[i][2].remove()
                        self.lines_super[i][2] = []
                    
                self.axes_record[order-1].grid(linestyle='dotted')
                self.axes_record[order-1].set_title('STFT: '+i)
                self.axes_record[order-1].set_xlabel('Time')
                self.axes_record[order-1].set_ylabel('Frequency')
            
            self.canvas.draw()
            self.sl_valuechange()
            
        def close_window():
            self.dialog.accept()    
            
        def ConceFT_cb_state(b):
                if b.isChecked() == True:
                    group_box4.setEnabled(True)
                else:
                    group_box4.setEnabled(False)
        
        
        
        
        #Confirm Data Loaded
        if not signal_dic:
            print ('data not loaded')
        
            
        #Create Dialog for settings
        self.dialog = QDialog()
        self.signal_list = QListWidget()
        for i in signal_dic:
            self.signal_list.addItem(i)
        vbox = QVBoxLayout(self.dialog)
        
        #Group Box1 for General Parameters
        group_box1 = QGroupBox()
        gb1 = QFormLayout()
        WindowLength = QLineEdit('377')
        WindowLength.setValidator(QIntValidator())
        WindowBandwidth = QLineEdit('10') 
        WindowBandwidth.setValidator(QIntValidator())
        tDS = QLineEdit('1')
        tDS.setValidator(QIntValidator())
        SuperImpose = QCheckBox()
        gb1.addRow("Window Length",WindowLength)
        gb1.addRow("Window Bandwidth",WindowBandwidth)
        gb1.addRow("tDS",tDS)
        gb1.addRow("Sumperimpose", SuperImpose)
        
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
        Close = QPushButton('Close')
        Close.clicked.connect(close_window)
        
        #Add Widgets and Set Layout
        vbox.addWidget(self.signal_list)
        vbox.addWidget(group_box1)
        if self.en_cb_2.currentIndex() != 0:
            vbox.addWidget(group_box2)
            vbox.addWidget(ConceFT_cb)
            vbox.addWidget(group_box4)
        if self.en_cb_2.currentIndex() == 1:
            vbox.addWidget(group_box3)
        vbox.addWidget(Run)
        vbox.addWidget(Close)
        
        #Exececute
        
        self.setLayout(vbox)
        self.dialog.setWindowTitle("STFT Settings")
        self.dialog.setWindowModality(Qt.ApplicationModal)
        self.dialog.exec_()
        
        
    def update_CWT(self):
        def plot_3():
            global signal_dic
            #Values for parameters
            LowFrequencyLimit = float(Low_freq.text())
            HighFrequencyLimit = float(High_freq.text())
            FAR = float(FrequencyAxisResolution.text())
            n_ConceFT = int(NoConceFT.text())
            ifSmooth = Smooth.isChecked()
            ifHemi = Hemi.isChecked()
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
            items = self.signal_list.selectedItems()
            for j in items:
                i = j.text(); order = signal_dic[i][6][3]
                self.if_superimpose[order-1] = SuperImpose.isChecked()
                SR = self.SR[i]; self.if_log[order-1] = self.zlog_cb3.currentIndex()
                signal_todo = signal_dic[i][3]
                t_todo = signal_dic[i][1]
                if self.en_cb_3.currentIndex() == 0:
                    z, _, _, _, y = ConceFT_CWT(t_todo, signal_todo, LowFrequencyLimit, HighFrequencyLimit, FAR, 0, opts_in, 0, 0) 
                elif self.en_cb_3.currentIndex() == 1:    
                    if ConceFT_cb.isChecked():
                        _, _, _, z, y = ConceFT_CWT(t_todo, signal_todo, LowFrequencyLimit, HighFrequencyLimit, FAR, n_ConceFT, opts_in, ifSmooth, ifHemi)  
                    else:
                        _, _, z, _, y = ConceFT_CWT(t_todo, signal_todo, LowFrequencyLimit, HighFrequencyLimit, FAR, 0, opts_in, 0, 0) 
            
                
                
                self.colorbars[order-1] = imageSQ(self.figure, self.axes_record[order-1], self.colorbars[order-1], t_todo, y*SR, np.abs(z), 99.5,self.if_log[order-1], 'Greys')
                self.axes_record[order-1].legend(shadow=True, fancybox=True)
        
            
                if self.if_superimpose[order-1]:
                    if not self.axes_super[order-1]:
                        self.axes_super[order-1] = self.axes_record[order-1].twinx()
                        self.axes_super[order-1].plot(signal_dic[i][1],signal_dic[i][3], '-', linewidth = 1,label=i,color='r')
                        self.lines_super[i][3] = self.axes_super[order-1].lines[-1]
                    else:
                        if self.lines_super[i][3]:
                            self.lines_super[i][3].remove()
                            self.lines_super[i][3] = []
                        self.axes_super[order-1].plot(signal_dic[i][1],signal_dic[i][3], '-', linewidth = 1,label=i,color='r')
                        self.lines_super[i][3] = self.axes_super[order-1].lines[-1]
                else:
                    if self.axes_super[order-1]and self.lines_super[i][3]:
                        self.lines_super[i][3].remove()
                        self.lines_super[i][3] = []
                    
                self.axes_record[order-1].grid(linestyle='dotted')
                self.axes_record[order-1].set_title('CWT: '+i)
                self.axes_record[order-1].set_xlabel('Time')
                self.axes_record[order-1].set_ylabel('Frequency')
            
            self.canvas.draw()
            self.sl_valuechange()
            
        def close_window():
            self.dialog.accept()  
            
        def ConceFT_cb_state(b):
                if b.isChecked() == True:
                    group_box4.setEnabled(True)
                else:
                    group_box4.setEnabled(False)
        
        
        class opts_in:
            motherwavelet = 'Cinfc'
            CENTER = 1
            FWHM = 0.3
        
        if not signal_dic:
            print ('data not loaded')
            return
        
        #Create Dialog for Settings
        self.dialog = QDialog()
        self.signal_list = QListWidget()
        for i in signal_dic:
            self.signal_list.addItem(i)
        vbox = QVBoxLayout(self.dialog)
        
        #Group Box1 for General Parameeters
        group_box1 = QGroupBox()
        
        gb1 = QFormLayout()
        SuperImpose = QCheckBox()
        gb1.addRow("Sumperimpose", SuperImpose)
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
        Close = QPushButton('Close')
        Close.clicked.connect(close_window)
        
        vbox.addWidget(self.signal_list)
        vbox.addWidget(MW_group_box)
        vbox.addWidget(group_box1)
        if self.en_cb_3.currentIndex() == 1:
            vbox.addWidget(group_box2)
            vbox.addWidget(ConceFT_cb)
            vbox.addWidget(group_box4)
        vbox.addWidget(Run)
        vbox.addWidget(Close)
        
        
        self.setLayout(vbox)
        self.dialog.setWindowTitle("CWT Settings")
        self.dialog.setWindowModality(Qt.ApplicationModal)
        self.dialog.exec_()
    
    def update_stran(self):
        def plot_4():
            global signal_dic
            #Values for parameters
            LowFrequencyLimit = float(Low_freq.text())
            HighFrequencyLimit = float(High_freq.text())
            FAR = float(FrequencyAxisResolution.text())
            ifSmooth = Smooth.isChecked()
            ifHemi = Hemi.isChecked()
            items = self.signal_list.selectedItems()
            for j in items:
                i = j.text(); order = signal_dic[i][6][4]
                self.if_superimpose[order-1] = SuperImpose.isChecked()
                SR = self.SR[i]; self.if_log[order-1] = self.zlog_cb4.currentIndex()
                signal_todo = signal_dic[i][3]
                t_todo = signal_dic[i][1] 
                if self.en_cb_4.currentIndex()==0:
                    z,y,_,_ = stran(signal_todo,t_todo,0)
                if self.en_cb_4.currentIndex()==1:
                    _,_,z,y = stran(signal_todo,t_todo,1,LowFrequencyLimit,HighFrequencyLimit,FAR,ifSmooth,ifHemi)
            
                self.colorbars[order-1] = imageSQ(self.figure, self.axes_record[order-1], self.colorbars[order-1], t_todo, y*SR, np.abs(z), 99.5,self.if_log[order-1], 'Greys')
                self.axes_record[order-1].legend(shadow=True, fancybox=True)
        
            
                if self.if_superimpose[order-1]:
                    if not self.axes_super[order-1]:
                        self.axes_super[order-1] = self.axes_record[order-1].twinx()
                        self.axes_super[order-1].plot(signal_dic[i][1],signal_dic[i][3], '-', linewidth = 1,label=i,color='r')
                        self.lines_super[i][4] = self.axes_super[order-1].lines[-1]
                    else:
                        if self.lines_super[i][4]:
                            self.lines_super[i][4].remove()
                            self.lines_super[i][4] = []
                        self.axes_super[order-1].plot(signal_dic[i][1],signal_dic[i][3], '-', linewidth = 1,label=i,color='r')
                        self.lines_super[i][4] = self.axes_super[order-1].lines[-1]
                else:
                    if self.axes_super[order-1] and self.lines_super[4][5]:
                        self.lines_super[i][4].remove()
                        self.lines_super[i][4] = []
                
                self.axes_record[order-1].grid(linestyle='dotted')
                self.axes_record[order-1].set_title('S Transform: '+i)
                self.axes_record[order-1].set_xlabel('Time')
                self.axes_record[order-1].set_ylabel('Frequency')
            
            self.canvas.draw()
            self.sl_valuechange()
            
        def close_window():
            self.dialog.accept()      
        
        if signal_dic == {}:
            print ('data not loaded')
            return
        
        #Create Dialog for Settings
        self.dialog = QDialog()
        self.signal_list = QListWidget()
        for i in signal_dic:
            self.signal_list.addItem(i)
        vbox = QVBoxLayout(self.dialog)
        
        #Group Box1 for General Parameeters
        group_box1 = QGroupBox()
        
        gb1 = QFormLayout()
        SuperImpose = QCheckBox()
        gb1.addRow("Sumperimpose", SuperImpose)
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
        Close = QPushButton('Close')
        Close.clicked.connect(close_window)
        
        vbox.addWidget(self.signal_list)
        vbox.addWidget(group_box1)
        if self.en_cb_4.currentIndex() == 1:
            vbox.addWidget(group_box2)
        vbox.addWidget(Run)
        vbox.addWidget(Close)
        
        self.setLayout(vbox)
        self.dialog.setWindowTitle("S Transform Settings")
        self.dialog.setWindowModality(Qt.ApplicationModal)
        self.dialog.exec_()
    
    def update_WVD(self):
        def plot_5():
            global signal_dic
            #Values for parameters
            FAR = float(FrequencyAxisResolution.text())
            LowFrequencyLimit = float(Low_freq.text())
            HighFrequencyLimit = float(High_freq.text())
            ifSmooth = Smooth.isChecked()
            ifHemi = Hemi.isChecked()
            items = self.signal_list.selectedItems()
            for j in items:
                i = j.text(); order = signal_dic[i][6][5]
                self.if_superimpose[order-1] = SuperImpose.isChecked()
                SR = self.SR[i]; self.if_log[order-1] = self.zlog_cb5.currentIndex()
                signal_todo = signal_dic[i][3]
                t_todo = signal_dic[i][1] 
                if self.en_cb_5.currentIndex()==0:
                    z,y,_,_ = wvd(signal_todo,None,0)
                if self.en_cb_5.currentIndex()==1:
                    _,_,z,y = wvd(signal_todo,None,1,FAR, 0, True, LowFrequencyLimit,HighFrequencyLimit,ifSmooth,ifHemi)
                
                self.colorbars[order-1] = imageSQ(self.figure, self.axes_record[order-1], self.colorbars[order-1], t_todo, y*SR, np.abs(z), 99.5,self.if_log[order-1], 'Greys')
                self.axes_record[order-1].legend(shadow=True, fancybox=True)
        
                if self.if_superimpose[order-1]:
                    if not self.axes_super[order-1]:
                        self.axes_super[order-1] = self.axes_record[order-1].twinx()
                        self.axes_super[order-1].plot(signal_dic[i][1],signal_dic[i][3], '-', linewidth = 1,label=i,color='r')
                        self.lines_super[i][5] = self.axes_super[order-1].lines[-1]
                    else:
                        if self.lines_super[i][5]:
                            self.lines_super[i][5].remove()
                            self.lines_super[i][5] = []
                        self.axes_super[order-1].plot(signal_dic[i][1],signal_dic[i][3], '-', linewidth = 1,label=i,color='r')
                        self.lines_super[i][5] = self.axes_super[order-1].lines[-1]
                else:
                    if self.axes_super[order-1]:
                        self.lines_super[i][5].remove() and self.lines_super[i][5]
                        self.lines_super[i][5] = []
           
                self.axes_record[order-1].grid(linestyle='dotted')
                self.axes_record[order-1].set_title('WVD: '+i)
                self.axes_record[order-1].set_xlabel('Time')
                self.axes_record[order-1].set_ylabel('Frequency')
            
            self.canvas.draw()
            self.sl_valuechange()
        
        
        def close_window():
            self.dialog.accept()      
        
        if signal_dic == {}:
            print ('data not loaded')
            return
        
        #Create Dialog for Settings
        self.dialog = QDialog()
        self.signal_list = QListWidget()
        for i in signal_dic:
            self.signal_list.addItem(i)
        vbox = QVBoxLayout(self.dialog)
        
        #Group Box1 for General Parameeters
        group_box1 = QGroupBox()
        gb1 = QFormLayout()
        SuperImpose = QCheckBox()
        gb1.addRow("Sumperimpose", SuperImpose)
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
        Close = QPushButton('Close')
        Close.clicked.connect(close_window)
        
        vbox.addWidget(self.signal_list)
        vbox.addWidget(group_box1)
        
        if self.en_cb_5.currentIndex() == 1:
            vbox.addWidget(group_box2)
        
        vbox.addWidget(Run)
        vbox.addWidget(Close)
        
        
        self.setLayout(vbox)
        self.dialog.setWindowTitle("WVD Settings")
        self.dialog.setWindowModality(Qt.ApplicationModal)
        self.dialog.exec_()
        
        
    def update_CFPH(self,mode):
        def plot_6():
            global signal_dic
            #Values for parameters
            class basicTF:
                win = int(WindowLength.text()) #4096
                hop = int(tDS.text()) #441;
                fs = int(Fs.text())
                fr = float(Fr.text());
                feat = Feat.currentText()
                    
            
            class advTF:
                num_tap = int(MT.text()) #num_tap(cc);
                win_type = WinType.currentText() # %{'Gaussian','Thomson','multipeak','SWCE'}; %}
                Smo = Smooth.isChecked()
                Rej = Rejection.isChecked()
                ths = float(Ths.text())
                HighFreq = float(High_freq.text())
                LowFreq = float(Low_freq.text())
                lpc = int(lpc_num.text())

            class cepR:
                g = float(G.text()) # % g(cc);%0.06;
                Tc = float(Tc_w.text())

            class P:
                num_s = int(Num_C.text())
                num_c = int(Num_S.text());
            
            if advTF.num_tap<1:
                print('MT must be larger than one')
                return
            
            items = self.signal_list.selectedItems()
            for j in items:
                i = j.text(); order = 1;self.SR[i] = basicTF.fs
                signal_todo = signal_dic[i][3]
                t_todo = signal_dic[i][1]
                z=[];y=[];plot_fn = 6
                
                if mode == "CEPS":
                    plot_fn = 6;order = signal_dic[i][6][6];self.if_log[order-1] = self.zlog_cb6.currentIndex()
                    _, z, _, _, _, _, y = CFPH(signal_todo, basicTF, advTF, cepR, P)
                elif mode == "TCEPS": 
                    plot_fn = 7; order = signal_dic[i][6][7];self.if_log[order-1] = self.zlog_cb7.currentIndex()
                    _, _, z, _, _, _, y = CFPH(signal_todo, basicTF, advTF, cepR, P)
                elif mode == "TFRR": 
                    plot_fn = 8; order = signal_dic[i][6][8];self.if_log[order-1] = self.zlog_cb8.currentIndex()
                    _, _, _, z, _, _, y = CFPH(signal_todo, basicTF, advTF, cepR, P)
                elif mode == "RTFR": 
                    plot_fn = 9; order = signal_dic[i][6][9]; self.if_log[order-1] = self.zlog_cb9.currentIndex()
                    _, _, _, _, z, _, y = CFPH(signal_todo, basicTF, advTF, cepR, P)
                
                self.if_superimpose[order-1] = SuperImpose.isChecked()
                
                
                
                    
                
                self.colorbars[order-1] = imageSQ(self.figure, self.axes_record[order-1], self.colorbars[order-1], t_todo, y*self.SR[i], np.abs(z), 99.5,self.if_log[order-1], 'Greys')
                self.axes_record[order-1].legend(shadow=True, fancybox=True)
        
                if self.if_superimpose[order-1]:
                    if not self.axes_super[order-1]:
                        self.axes_super[order-1] = self.axes_record[order-1].twinx()
                        self.axes_super[order-1].plot(signal_dic[i][1],signal_dic[i][3], '-', linewidth = 1,label=i,color='r')
                        self.lines_super[i][plot_fn] = self.axes_super[order-1].lines[-1]
                    else:
                        if self.lines_super[i][plot_fn]:
                            self.lines_super[i][plot_fn].remove()
                            self.lines_super[i][plot_fn] = []
                        self.axes_super[order-1].plot(signal_dic[i][1],signal_dic[i][3], '-', linewidth = 1,label=i,color='r')
                        self.lines_super[i][plot_fn] = self.axes_super[order-1].lines[-1]
                else:
                    if self.axes_super[order-1] and self.lines_super[i][plot_fn]:
                        self.lines_super[i][plot_fn].remove()
                        self.lines_super[i][plot_fn] = []
           
                self.axes_record[order-1].grid(linestyle='dotted')
                self.axes_record[order-1].set_xlabel('Time')
                
                if plot_fn!= 6:self.axes_record[order-1].set_ylabel('Frequency')
                else:self.axes_record[order-1].set_ylabel('quefrency')
                
                if plot_fn == 6:self.axes_record[order-1].set_title('Ceps: '+i)
                elif plot_fn == 7:self.axes_record[order-1].set_title('tCeps: '+i)
                elif plot_fn == 8:self.axes_record[order-1].set_title('TFRR: '+i)
                elif plot_fn == 9:self.axes_record[order-1].set_title('RTFR: '+i)
            
            self.canvas.draw()
            self.sl_valuechange()
        
        
        def close_window():
            self.dialog.accept()      
        
        
            
        if signal_dic == {}:
            print ('data not loaded')
            return
        
        #Create Dialog for Settings
        self.dialog = QDialog()
        self.signal_list = QListWidget()
        for i in signal_dic:
            self.signal_list.addItem(i)
        vbox = QVBoxLayout(self.dialog)
        
        #Group Box1 for General Parameeters
        group_box1 = QGroupBox()
        gb1 = QFormLayout()
        WindowLength = QLineEdit('377')
        WindowLength.setValidator(QIntValidator())
        WinType = QComboBox()
        WinType.addItems(['Gauss','Hamming','RECTANG','RECT','HANNING','HANN','KAISER','NUTTALL','BLACKMAN','BLACKMANHARRIS','HARRIS','BARTLETT','TRIANG','BARTHANN','PAPOULIS','PARZEN','HANNA','DOLPH','DOLF','NUTBESS','SPLINE','FLATTOP','POISSON','HANNINGPOISSON','CAUCHY'])
        tDS = QLineEdit('1')
        tDS.setValidator(QIntValidator())
        Fs = QLineEdit('50')
        Fs.setValidator(QIntValidator())
        Fr = QLineEdit('0.02')
        Fr.setValidator(QDoubleValidator())
        Feat = QComboBox()
        Feat.addItems(["STFT","SST11"])
        MT = QLineEdit('1'); MT.setValidator(QIntValidator())
        Smooth = QCheckBox()
        Rejection = QCheckBox()
        Ths = QLineEdit('1E-9') ; Ths.setValidator(QDoubleValidator())
        High_freq = QLineEdit(str(8/50)); High_freq.setValidator(QDoubleValidator())
        Low_freq = QLineEdit(str(0.1/50));Low_freq.setValidator(QDoubleValidator())
        lpc_num = QLineEdit('0');lpc_num.setValidator(QIntValidator())
        SuperImpose = QCheckBox()
        
        gb1.addRow("Window Length",WindowLength)
        gb1.addRow("Window Type", WinType)
        gb1.addRow("tDS",tDS)
        gb1.addRow("fs", Fs)
        gb1.addRow("fr", Fr)
        if mode == "TFRR":
            gb1.addRow("Feat",Feat)
        if mode in ["TFRR","TFR","RTFR"]:
            gb1.addRow("MT", MT)
            gb1.addRow("Smooth", Smooth)
            gb1.addRow("Rej", Rejection)
            gb1.addRow("Ths", Ths)
        
        gb1.addRow("High Frequency", High_freq)
        gb1.addRow("Low Frequency", Low_freq)
        gb1.addRow("lpc_num", lpc_num)
        gb1.addRow("Sumperimpose", SuperImpose)
        group_box1.setLayout(gb1)
        
        #sq box
        group_box2 = QGroupBox()
        gb2 = QFormLayout()
        G = QLineEdit('0.1')
        G.setValidator(QDoubleValidator())
        Tc_w = QLineEdit('0')
        Tc_w.setValidator(QDoubleValidator())
        Num_C = QLineEdit('1'); Num_C.setValidator(QIntValidator())
        Num_S = QLineEdit('1'); Num_S.setValidator(QIntValidator())
        gb2.addRow("G", G)
        gb2.addRow("Tc", Tc_w)
        gb2.addRow("num_c", Num_C)
        if mode in ["TFRR","TFR","RTFR"]:
            gb2.addRow("num_c", Num_S)
        group_box2.setLayout(gb2)
        
        #Run Button
        Run = QPushButton('Run')
        Run.clicked.connect(plot_6)
        Close = QPushButton('Close')
        Close.clicked.connect(close_window)
        
        vbox.addWidget(self.signal_list)
        vbox.addWidget(group_box1)
        vbox.addWidget(group_box2)
        
        vbox.addWidget(Run)
        vbox.addWidget(Close)
        
        
        self.setLayout(vbox)
        self.dialog.setWindowTitle("Ceps Settings")
        self.dialog.setWindowModality(Qt.ApplicationModal)
        self.dialog.exec_()
        
def update_time_lim():
    global t_min, t_max
    if t_min_dic and t_max_dic:
        t_min = min(t_min_dic.values())
        t_max = max(t_max_dic.values())
    print(t_min,t_max)

def main():
   app = QApplication(sys.argv)
   ex = Window()
   ex.showMaximized()
   sys.exit(app.exec_())
	
if __name__ == '__main__':
   main()       
        