# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 09:41:37 2017
This program generates plots for 1D SWE program.

@author: VL Srinivas (https://www.linkedin.com/in/vl-srinivas/)
"""

from numpy import loadtxt
#from pylab import figure, ioff
from matplotlib import pyplot as pp
#from matplotlib import figure
from glob import glob
import os

#ioff()

fileNames = glob("/home/srinivas/Desktop/vls/files/h_*.dat")
fileNames.sort()

for fn in fileNames:
    #fig = figure
    #pyplot.figure(figsize=(10,6))
    pp.figure(figsize=(10,8))   
    f = file(fn, 'r')
    nx = int(f.readline().split(":")[1])
    #nt = int(f.readline().split(":")[2])
    f.close()
    
    x,h = loadtxt(fn, skiprows=2, unpack=True)
    pp.plot(x, h, lw = 3)
    pp.xlim(0.0, 1.0)
    pp.ylim(-2.0, 5.0)
    pp.title("1D Shallow water equation")
    pp.savefig(fn.replace(".dat",".png"))
    print fn
#os.system("eog " + fileNames[0].replace(".dat",".png"))
    
    
