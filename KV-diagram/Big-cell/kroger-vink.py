# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 21:22:28 2017

@author: minhdinh
"""


import sys
from PyQt5 import QtCore, QtGui, uic, QtWidgets
from decimal import Decimal
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)

import numpy
#import os, sys
import csv
import scipy.optimize 
#from franges import frange
#from sympy.solvers.solvers import nsolve
#from sympy import Symbol
#import plotly.plotly as pltly      # plotting functions
#import plotly.tools as tls         # plotly tools
#import plotly.graph_objs as go     # plot and configuration tools : Scatter, Line, Layout


qtCreatorFile = "K-V.ui" # Enter file here.

Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

class Ui(QtWidgets.QMainWindow):
    def __init__(self):
        super(Ui, self).__init__() # Call the inherited classes __init__ method
        uic.loadUi(qtCreatorFile, self) # Load the .ui file
        
        #self.button = self.findChild(QtWidgets.QPushButton, 'calculateButton') # Find the button
        self.calculateButton.clicked.connect(self.calculateButtonPressed) # Remember to pass the definition/method, not the return value
        self.DefectFileFormatInput.setText("defects-*.txt")
        self.doscarLocation.setText("./DOSCAR")
        self.dopingConcentrationInput.setText("10000")
        self.minOxygenP.setText("1e-20")
        self.maxOxygenP.setText("1")
        #self.defectListInput = self.findChild(QtWidgets.QLineEdit, 'DefectListInput')
        self.actionExit.triggered.connect(self.close)
        self.temperatureSlider.valueChanged.connect(self.valuechange)
        
        
        

        
        
        # CONSTANT
        self.k = 8.6173303*10**-5 #eV/K
        self.temperature = 600 # K
        #mu_e = 0.92566997 # eV
        self.dielectric_const = 24;
        #E_MK =  2.337359 # eV
        self.E_MK = 1.907001 # large cell
        
        # Tabulated info in thermochemical table: 
        # Oxygen chemical potential at p0=1atm at different T (K)
        #   PHYSICAL REVIEW B 65 035406
        self.mu_H2O_p0_table = [[298.15, -0.48086], [300, -0.48448], [400, -0.68568],\
        [500, -0.89585], [600, -1.1134], [700, -1.3372], [800, -1.5666],\
        [900, -1.801], [1000, -2.04], [1100, -2.2833], [1200, -2.5306],\
        [1300, -2.7817], [1400, -3.0364], [1500, -3.2945], [1600, -3.5558],\
        [1700, -3.8203], [1800, -4.0877], [1900, -4.358], [2000, -4.6311]]
        self.mu_O_p0_table = [[298.150000, -0.271962], [300.000000, -0.273931],\
        [350.000000, -0.327734], [400.000000, -0.382644],\
        [450.000000, -0.438528], [500.000000, -0.495295],\
        [600.000000, -0.611189], [700.000000, -0.729865],\
        [800.000000, -0.850985], [900.000000, -0.974297],\
        [1000.000000, -1.09959], [1100.000000, -1.22669],\
        [1200.000000, -1.35546], [1300.000000, -1.48577],\
        [1400.000000, -1.61752], [1500.000000, -1.7506],\
        [1600.000000, -1.88495], [1700.000000, -2.02049],\
        [1800.000000, -2.15717], [1900.000000, -2.29491],\
        [2000.000000, -2.43368]]
        
        
        #mu_O_p0_table = [[100,-0.08],[200,-0.17],[300,-0.27],[400,-0.38],[500,-0.50],\
        #           [600,-0.61],[700,-0.73],[800,-0.85],[900,-0.98],[1000,-1.10]]
        
        # Energy constant:
        # MxOy
        self.xM = 1
        self.yO = 2
        
        self.E_Metal = -3.8162 # eV. Bulk metal
        self.E_MetalOxide = -19.6874 # eV. Energy per unit cell
        self.E_water = -14.2176 # eV. One water molecule
        
        E_O2_over = 1.22369
        self.E_O2 = -9.8579 + E_O2_over # eV. Energy of O2 molecule
        #E_MetalFluoride = -23.854444 # eV 
        self.dopingConcentration = 0.0001 
        
        self.cellsize = 32 # formula units
        self.dosCellsize = 16
        
        
        

        
        
        
        
        
        
        self.show() # Show the GUI
        
    def valuechange(self):
      size = self.temperatureSlider.value()
      size = (size - 1)*100 + 300
      self.temperature = size
      self.temperatureInput.setText(str(size)+" K")
        
    def calculateButtonPressed(self):
        # This is executed when the button is pressed
        
        # READ DOSCAR
        outputText =  "Calculating ..."
        self.Output.setText(outputText)
        doscarFile = self.doscarLocation.text()
        file = open(doscarFile,"r")
        reader = csv.reader(file,delimiter=' ')
        
        doscar = []
        for row in reader:       
            doscar.append([x for x in row if x])
        #header = doscar[5]
        self.doscarData = numpy.array(doscar[6:(int(doscar[5][2])+6)]).astype(numpy.float)
        #gap,vbm,cbm = get_bandgap(doscarFile,1e-3)
        self.gap, self.vbm, self.cbm = [3.6, 6.2704, 3.6+6.2704]
        file.close()
        #print(self.gap,self.vbm,self.cbm)

        
        defectFormat = self.DefectFileFormatInput.text()
        defectName = self.defectListInput.text()
        defectFile = defectFormat.replace("*",defectName)
        file = open(defectFile,"r")
        reader = csv.reader(file,delimiter=' ')
        print(defectFile)
        self.defect = []
        for row in reader:       
            self.defect.append([x for x in row if x])
        file.close()
        self.totalDopingDefects = numpy.count_nonzero(numpy.char.find(numpy.array(self.defect)[:,4],'1')+1)
        
        #self.temperature = int(self.TemperatureInput.value())
        
        self.mu_O_p0 = [x[1] for x in self.mu_O_p0_table if x[0] == self.temperature][0]
        self.mu_H2O_p0 = [x[1] for x in self.mu_H2O_p0_table if x[0] == self.temperature][0]
        self.mu_O_min = (self.E_MetalOxide - self.xM*self.E_Metal) / self.yO
        self.mu_O_max = (self.E_O2) / 2 + self.mu_O_p0
        self.mu_H = 1/2*(self.E_water + self.mu_H2O_p0 - 1/2*self.E_O2 - self.mu_O_p0)
        #minOP = float(self.minOxygenP.getText())
        #maxOP = float(self.maxOxygenP.getText())
        output = []
        for j in numpy.arange(self.mu_O_max, self.mu_O_min, -0.2):         
            self.mu_O = j
            pO2 = numpy.exp(2*(self.mu_O-self.E_O2/2-self.mu_O_p0)/self.k/self.temperature)
            self.mu_Metal = 1/self.xM*(self.E_MetalOxide - self.yO*self.mu_O)         
            #self.Output.setDisabled(True)
            outputText += "\n pO2 = "+ "{:.2E}".format(Decimal(pO2))
            self.Output.setText(outputText)
            print(str(pO2))
            [mu_e2, mu_doping2] = scipy.optimize.fsolve(self.dopingLevelandTotalCharge, (2,-10), maxfev = 20000)
            temp_nc = self.nc(self.doscarData, self.temperature, mu_e2+self.vbm, self.cbm)
            temp_pv = self.pv(self.doscarData, self.temperature, mu_e2+self.vbm, self.vbm)
            tempOutput = [pO2, mu_e2, temp_nc, temp_pv]
            for i in range(1,len(self.defect)):
                charge = int(self.defect[i][0])
                nM = [int(x) for x in self.defect[i][1:5]]
            
                noSite = int(self.defect[i][-1])
            
                E = float(self.defect[i][-2])
                mu = [self.mu_O, self.mu_Metal, self.mu_H, mu_doping2]
                #print([E, charge, vbm, mu_e, nM, mu])
                E_form = self.formation_energy (E, charge, self.vbm, mu_e2, nM, mu)
                c = self.concentration(E_form, self.temperature, noSite)
                tempOutput.append(c)
            output.append(tempOutput)
            
        
            #print(mu_e2)
            #self.Output.setDisabled(False)
            #self.Output.setText("{:.2f}".format(mu_e2))   


    def concentration(self, Eform, temperature, noSite):
        return noSite*numpy.exp(-Eform/(self.k*temperature))

    def formation_energy (self, Etot, charge, vbm, mu_e, nMu, mu): 
        energy = Etot - self.E_MetalOxide * self.cellsize \
        + charge**2*self.E_MK/self.dielectric_const + charge*(vbm+mu_e)  \
        - sum([x*y for x,y in zip(nMu,mu)])
        #print([Etot, charge, vbm, mu_e, nMu, mu])
        #print(energy)
        #wait = input("PRESS ENTER TO CONTINUE.")
        return energy
    
    
    
    def totalDefectsCharge (self, mu_e,temperature, vbm, mu):
        defect = self.defect
        tot_charge = 0 
        for i in range(1,len(defect)):
            charge = int(defect[i][0])
            nM = [int(x) for x in defect[i][1:5]]
            
            noSite = int(defect[i][-1])
            
            E = float(defect[i][-2])
            #print([E, charge, vbm, mu_e, nM, mu])
            E_form = self.formation_energy (E, charge, vbm, mu_e, nM, mu)
            c = self.concentration(E_form, temperature, noSite)
            tot_charge += charge*c            
        return tot_charge
    
    
    def dopingLevel(self, variables):
        defect = self.defect
        (mu_e, mu_Doping) = variables
        dopingC = 0
        mu = [self.mu_O, self.mu_Metal, self.mu_H, mu_Doping]
        for i in range(-self.totalDopingDefects,0):
            charge = int(defect[i][0])
            nM = [int(x) for x in defect[i][1:5]]
            
            noSite = int(defect[i][-1])
            
            E = float(defect[i][-2])
            E_form = self.formation_energy (E, charge, self.vbm, mu_e, nM, mu)
            #print(E_form)       
            c = self.concentration(E_form, self.temperature, noSite)
            dopingC += c
        return c
    
    def dopingLevelandTotalCharge (self, variables):
        defect = self.defect
        (mu_e, mu_Doping) = variables
        dopingC = 0
        mu = [self.mu_O, self.mu_Metal, self.mu_H, mu_Doping]
        for i in range(-self.totalDopingDefects,0):
            charge = int(defect[i][0])
            nM = [int(x) for x in defect[i][1:5]]
            
            noSite = int(defect[i][-1])
            
            E = float(defect[i][-2])
            E_form = self.formation_energy (E, charge, self.vbm, mu_e, nM, mu)
            #print(E_form)       
            c = self.concentration(E_form, self.temperature, noSite)
            dopingC += c
        #print(dopingC)    
        firstEQ = dopingC - self.dopingConcentration
        
        secondEQ = -self.nc(self.doscarData, self.temperature, mu_e+self.vbm, self.cbm) \
            + self.pv(self.doscarData, self.temperature, mu_e+self.vbm, self.vbm) \
            + self.totalDefectsCharge (mu_e,self.temperature, self.vbm, mu)
        #print([mu_e, mu_Doping, firstEQ, secondEQ])
        return numpy.asarray([firstEQ, secondEQ])
        
    
    def get_bandgap(location,tol):
        doscar = open(location)
        for i in range(6):
            l=doscar.readline()
        efermi = float(l.split()[3])
        step1 = doscar.readline().split()[0]
        step2 = doscar.readline().split()[0]
        step_size = float(step2)-float(step1)
        not_found = True
        while not_found:
            l = doscar.readline().split()
            e = float(l.pop(0))
            dens = 0
            for i in range(int(len(l)/2)):
                dens += float(l[i])
            if e < efermi and dens > tol:
                bot = e
            elif e > efermi and dens > tol:
                top = e
                not_found = False
        if top - bot < step_size*2:
            return 0,0,0
        else:
            #return top - bot,bot-efermi,top-efermi
            return top - bot,bot,top
        
    def fermi_dirac(self, E, temperature, mu_e ):
    
        return 1/(1+numpy.exp((E-mu_e)/(self.k*temperature)))
    
    def nc (self, doscar, temperature, fermi, cbm ):
       "Conduction band electron"
       temp = [x for x in doscar if x[0]>=cbm]
       #step = (temp[-1][0]-temp[0][0])/len(temp) 
       sum = 0
       prevNumberOfE = temp[0][2]
       for row in temp:
           sum += (row[2]-prevNumberOfE)*self.fermi_dirac (row[0],temperature,fermi)
           prevNumberOfE = row[2]
    #   print(sum)
       return sum / self.dosCellsize
    
    def pv (self, doscar, temperature, fermi, vbm ):
        "Valence band hole"
        temp = [x for x in doscar if x[0]<=vbm]
        #step = (temp[-1][0]-temp[0][0])/len(temp)   
        sum = 0
        prevNumberOfE = temp[0][2]
        for row in temp:
            sum += (row[2]-prevNumberOfE)*self.fermi_dirac(-row[0],temperature,-fermi)
            prevNumberOfE = row[2]
    #    print(sum)
        return sum / self.dosCellsize
   
    def addmpl(self, fig):
        self.canvas = FigureCanvas(fig)
        self.mplvl.addWidget(self.canvas)
        self.canvas.draw()
#def constrainedFunction(x, f, lower, upper, minIncr=0.001):
#     x = numpy.asarray(x)
#     lower = numpy.asarray(lower)
#     upper = numpy.asarray(upper)
#     xBorder = numpy.where(x<lower, lower, x)
#     xBorder = numpy.where(x>upper, upper, xBorder)
#     fBorder = f(xBorder)
#     #print(fBorder)
#     distFromBorder = (numpy.sum(numpy.where(x<lower, lower-x, 0.))
#                      +numpy.sum(numpy.where(x>upper, x-upper, 0.)))
#     return (fBorder + (fBorder+numpy.where(fBorder>0, minIncr, -minIncr))*distFromBorder)   


if __name__ == "__main__":
    numpy.seterr(all='ignore')
    


#    nc(doscarData, temperature, mu_e+vbm, cbm)
#    pv(doscarData, temperature, mu_e+vbm, vbm)
    
    
    ### DEFECTS CALCULATION:
    
    # READ defect list
    # defectFile = "defects.txt"
    #defectFile = "defects-Ti.txt"
    #file = open(defectFile,"r")
    #reader = csv.reader(file,delimiter=' ')
    
    #defect = []
    #for row in reader:       
    #    defect.append([x for x in row if x])
    
    #totalDopingDefects = numpy.count_nonzero(numpy.char.find(numpy.array(defect)[:,4],'1')+1)


    
    #mu_H = (E_water - mu_O)/2
    #pO2 = numpy.exp(2*(mu_O-E_O2/2-mu_O_p0)/k/temperature)
    #print(pO2)
    
    
    "SOLVING FOR MU_E"
    #[mu_e2, mu_Doping2] = scipy.optimize.fsolve(dopingLevelandTotalCharge, (2,-10), maxfev = 20000)
    #print(mu_e2)
    
    
    #[mu_e, mu_Doping] = scipy.optimize.broyden1(dopingLevelandTotalCharge,(1,-10),f_tol=1e-14)    
    #[mu_e2, mu_Doping2] = scipy.optimize.root(constrainedFunction, numpy.asarray([1,-10]),  args=(dopingLevelandTotalCharge, numpy.asarray([0., -30.]), numpy.asarray([3., 0.])))
    #[mu_e, mu_Doping] = scipy.optimize.broyden1(dopingLevelandTotalCharge, (1.4,-6), f_tol = 1e-10)
    #[mu_e2, mu_Doping2] = scipy.optimize.root(dopingLevelandTotalCharge, [1,-10], method = 'hybr')
    app = QtWidgets.QApplication(sys.argv)
    window = Ui()
    
    window.show()
    sys.exit(app.exec_())
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    