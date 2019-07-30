##slices should be there already, named as S0.01,S0.02, etc
import numpy as np
import os
#import sys
#sys.modules[__name__].__dict__.clear()
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
##############################
##############################
#############
Smin=0.01
Smax=0.6
deltaS=0.01
nS=round((Smax-Smin)/deltaS+1)
print nS
SS=np.linspace(Smin,Smax,nS)
print SS
print len(SS)
L=0.4976


####### initializing files #######
path1='I:/PostProcess/pvout/'
path2=".csv"



index=1
for dS in SS:
    
       
    Sa =  St + dS*sh  #absolute S
    path = path1 + caseName + "_S" + str(np.around(dS)) + path2
    print "Current Location + " + str(Sa)
    source = FindSource("S" + str(np.around(dS,2)))
    
    # save data
    SaveData(path, proxy=source)
    # get active source.
    
    