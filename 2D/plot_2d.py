# -*- coding: utf-8 -*-
"""
MIT License
Copyright (c) 2017 VL Srinivas

Created on Fri Sep 29 09:41:37 2017
This program generates plots for 2D SWE program.

@author: VL Srinivas (https://www.linkedin.com/in/vl-srinivas/)
"""
import numpy as np
from numpy import loadtxt
#from pylab import figure, ioff
from matplotlib import pyplot as pp
from matplotlib import rcParams
#from matplotlib import figure
from glob import glob
import os
import time
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 12

#ioff()
#start = time.time()

fileNames = glob("/home/srinivas/Desktop/vls_2/files/h_*.dat")
fileNames.sort()

for fn in fileNames:
    #fig = figure
    pp.figure(figsize=(10,6))
#    pp.figure(figsize=(10,8))   
    f = file(fn, 'r')
    nx = int(f.readline().split(":")[1])
    ny = int(f.readline().split(":")[1])    
    nt = int(f.readline().split(":")[1])
    time = float(f.readline().split(":")[1])
    f.close()
    
    # Loading the text and reshaping
    x,y,h = loadtxt(fn, skiprows= 8, unpack=True)
    x = np.reshape(x,(nx,ny))
    y = np.reshape(y,(nx,ny))
    h = np.reshape(h,(nx,ny))
    pp.cla()
    pp.clf()
    pp.xticks([])
    pp.yticks([])
    
    pp.contourf(x, y, h, 100)
    pp.axes().set_aspect("equal")
    pp.xlabel('x')
    pp.ylabel('y')
#    pp.xlim(0.0, 1.0)
#    pp.ylim(-2.0, 5.0)
    pp.title("2D Shallow water equation, Time=%5.3f"%time)
    pp.colorbar()
    pp.savefig(fn.replace(".dat",".png"))
    print fn
#os.system("eog " + fileNames[0].replace(".dat",".png"))

#elapsed = float( time.time() - start)
#print('Time taken:', elapsed,'sec') 
    
