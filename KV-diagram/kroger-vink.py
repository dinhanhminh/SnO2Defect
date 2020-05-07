# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 21:22:28 2017

@author: minhdinh
"""

import numpy
#import os, sys
import csv
import scipy.optimize 
from franges import frange
#import plotly.plotly as pltly      # plotting functions
#import plotly.tools as tls         # plotly tools
#import plotly.graph_objs as go     # plot and configuration tools : Scatter, Line, Layout

def concentration( Eform, temperature, noSite):
    return noSite*numpy.exp(-Eform/(k*temperature))

def formation_energy (Etot, charge, vbm, mu_e, nMu, mu): 
    energy = Etot - E_MetalOxide*cellsize \
    + charge**2*E_MK/dielectric_const + charge*(vbm+mu_e)  \
    - sum([x*y for x,y in zip(nMu,mu)])
    #print([Etot, charge, vbm, mu_e, nMu, mu])
    #print(energy)
    #wait = input("PRESS ENTER TO CONTINUE.")
    return energy



def totalDefectsCharge (mu_e,temperature, vbm, mu):
    tot_charge = 0 
    for i in range(1,len(defect)):
        charge = int(defect[i][0])
        nM = [int(x) for x in defect[i][1:5]]
        
        noSite = int(defect[i][-1])
        
        E = float(defect[i][-2])
        #print([E, charge, vbm, mu_e, nM, mu])
        E_form = formation_energy (E, charge, vbm, mu_e, nM, mu)
        c = concentration(E_form, temperature, noSite)
        tot_charge += charge*c
        
        
        
    return tot_charge
def dopingLevel(variables):
    (mu_e, mu_Doping) = variables
    dopingC = 0
    mu = [mu_O, mu_Metal, mu_H, mu_Doping]
    for i in range(-totalDopingDefects,0):
        charge = int(defect[i][0])
        nM = [int(x) for x in defect[i][1:5]]
        
        noSite = int(defect[i][-1])
        
        E = float(defect[i][-2])
        E_form = formation_energy (E, charge, vbm, mu_e, nM, mu)
        #print(E_form)       
        c = concentration(E_form, temperature, noSite)
        dopingC += c
    return c

def dopingLevelandTotalCharge (variables):
    (mu_e, mu_Doping) = variables
    
    dopingC = 0
    mu = [mu_O, mu_Metal, mu_H, mu_Doping]
    for i in range(-totalDopingDefects,0):
        charge = int(defect[i][0])
        nM = [int(x) for x in defect[i][1:5]]
        
        noSite = int(defect[i][-1])
        
        E = float(defect[i][-2])
        E_form = formation_energy (E, charge, vbm, mu_e, nM, mu)
        #print(E_form)       
        c = concentration(E_form, temperature, noSite)
        dopingC += c
    #print(dopingC)    
    firstEQ = dopingC - dopingConcentration
    
    secondEQ = -nc(doscarData, temperature, mu_e+vbm, cbm) \
        + pv(doscarData, temperature, mu_e+vbm, vbm) \
        + totalDefectsCharge (mu_e,temperature, vbm, mu)
    #print([firstEQ, secondEQ])
    return [firstEQ, secondEQ]
    

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
    
def fermi_dirac( E, temperature, mu_e ):

    return 1/(1+numpy.exp((E-mu_e)/(k*temperature)))

def nc ( doscar, temperature, fermi, cbm ):
   "Conduction band electron"
   temp = [x for x in doscar if x[0]>=cbm]
   #step = (temp[-1][0]-temp[0][0])/len(temp) 
   sum = 0
   prevNumberOfE = temp[0][2]
   for row in temp:
       sum += (row[2]-prevNumberOfE)*fermi_dirac (row[0],temperature,fermi)
       prevNumberOfE = row[2]
#   print(sum)
   return sum

def pv ( doscar, temperature, fermi, vbm ):
    "Valence band hole"
    temp = [x for x in doscar if x[0]<=vbm]
    #step = (temp[-1][0]-temp[0][0])/len(temp)   
    sum = 0
    prevNumberOfE = temp[0][2]
    for row in temp:
        sum += (row[2]-prevNumberOfE)*fermi_dirac(-row[0],temperature,-fermi)
        prevNumberOfE = row[2]
#    print(sum)
    return sum
   
    


if __name__ == "__main__":
    numpy.seterr(all='ignore')
    
    # READ DOSCAR
    
    doscarFile = "DOSCAR"
    file = open(doscarFile,"r")
    reader = csv.reader(file,delimiter=' ')
    
    doscar = []
    for row in reader:       
        doscar.append([x for x in row if x])
    header = doscar[5]
    doscarData = numpy.array(doscar[6:(int(doscar[5][2])+6)]).astype(numpy.float)
    #gap,vbm,cbm = get_bandgap(doscarFile,1e-3)
    gap, vbm, cbm = [3.6, 6.2704, 3.6+6.2704]
    file.close()
#    print(gap,vbm,cbm)
    
    
    # READ defect list
#    defectFile = "defects.txt"
    defectFile = "defects-Zr.txt"
    file = open(defectFile,"r")
    reader = csv.reader(file,delimiter=' ')
    
    defect = []
    for row in reader:       
        defect.append([x for x in row if x])
    
    
    # CONSTANT
    k = 8.6173303*10**-5 #eV/K
    temperature = 300 # K
    mu_e = 0.92566997 # eV
    dielectric_const = 24;
    E_MK =  2.337359 # eV
    
    # Tabulated info in thermochemical table: 
    # Oxygen chemical potential at p0=1atm at different T (K)
    #   PHYSICAL REVIEW B 65 035406
    mu_O_p0_table = [[100,-0.08],[200,-0.17],[300,-0.27],[400,-0.38],[500,-0.50],\
               [600,-0.61],[700,-0.73],[800,-0.85],[900,-0.98],[1000,-1.10]]
    mu_O_p0 = [x[1] for x in mu_O_p0_table if x[0] == temperature][0]
    # Energy constant:
    # MxOy
    xM = 1
    yO = 2
    
    E_Metal = -3.8162 # eV. Bulk metal
    E_MetalOxide = -19.6874 # eV. Energy per unit cell
    E_water = -14.2176 # eV. One water molecule
    E_O2 = -9.8579 # eV. Energy of O2 molecule
    E_MetalFluoride = -23.854444 # eV 
    dopingConcentration = 0.0001 
    totalDopingDefects = numpy.count_nonzero(numpy.char.find(numpy.array(defect)[:,4],'1')+1)
    
    mu_O_min = (E_MetalOxide - xM*E_Metal) / yO
    mu_O_max = (E_O2) / 2
    
    mu_O = mu_O_p0
#    nc(doscarData, temperature, mu_e+vbm, cbm)
#    pv(doscarData, temperature, mu_e+vbm, vbm)
    
    cellsize = 16 # formula units
    ### DEFECTS CALCULATION:
    
    
    #for mu_O in frange(mu_O_min, mu_O_max, 1):
    mu_O = -5.5
    mu_Metal = 1/xM*(E_MetalOxide - yO*mu_O)
    mu_H = (E_water - mu_O)/2
    
    [mu_e, mu_Doping] = scipy.optimize.fsolve(dopingLevelandTotalCharge,(0,-1))
    pO2 = numpy.exp((mu_O-E_O2/2-mu_O_p0)/k/temperature)
    print(pO2, mu_e)
     

#        mu_F = (E_MetalFluoride - mu_Metal)/4 # SnF4 
#        mu = [mu_O, mu_Metal, mu_H, mu_Doping]
        # NET NEUTRALITY EQUATION:
#        net_neutral = lambda x : -nc(doscarData, temperature, x+vbm, cbm) \
#        + pv(doscarData, temperature, x+vbm, vbm) \
#        + totalDefectsCharge (x,temperature, vbm)
#        mu_e = fsolve(net_neutral, 1)
#        print(mu_O,mu_e[0])        
#    totalDefectsCharge(1, 300,vbm)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    