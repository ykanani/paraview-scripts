import numpy as np
import os
import struct
import array
#import sys
#sys.modules[__name__].__dict__.clear()
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
##############################
##############################
L=0.4976

###############################################
n= 50
delta=0.02
zmin=0.1#0.117 #0.1
zmax=0.2#0.137 #0.154
nZ=2
Z=np.linspace(Smin,Smax,nZ)
print Z
Pr = 0.711
nu = 1.55/100000.0
alfa = nu/Pr
gradT = 38000.0
T0=20
Uexit = 15.93
Uinlet = 4
#############
###############################################
z=0

###############################
## find source
sourceName = 'ps10' #raw_input("Please enter source name: ")
var='T1'
caseName = sourceName # raw_input("Please enter case name: ")
print sourceName
source = FindSource(sourceName + '.foam')
# get active view
renderView1 = GetActiveViewOrCreate('RenderView')



#####################################surface

# create a new 'Slice' on the surface

slice1 = Slice(Input=source)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [0, 0, 0.01]
slice1.SliceType.Normal = [0,0, 1]


##############################
##############################

####### initializing files #######
XX=[]
YY=[]
ZZ=[]
for z in Z:
    
    
    print "Current Location + " + str(z)
    
    # modify slice 2
    slice1.SliceType.Origin = [0, 0, z]
    slice1.SliceOffsetValues = [0.0, 0.01]
    slice1.SliceType.Normal = [0,0, 1]
    
    pdi = servermanager.Fetch(slice1)
    print pdi
    print pdi.GetBlock(0)
    Npoints = pdi.GetBlock(0).GetNumberOfPoints()
    print Npoints
    for i in range(Npoints):
        XX.append(pdi.GetBlock(0).GetPoint(1)[0])
        YY.append(pdi.GetBlock(0).GetPoint(1)[1])
        ZZ.append(pdi.GetBlock(0).GetPoint(1)[2])
    VarArray = pdi.GetBlock(0).GetPointData().GetArray(var)
    print VarArray
    s = ''
    s = s.join([struct.pack('f',VarArray.GetValue(i)) for i in range(T1Array.GetNumberOfTuples())])

    # now make array
    a = array.array('f')
    a.fromstring(s)
    print np.asarray(a) 
    
# check
    #print len(a), a
    m=T1Array.GetNumberOfTuples()
    T1=[]
    alaki
    #print T1Mean.GetTuple(0)[0]
    # for j in range(0,m):
        # T.append(T1Mean.GetValue(j))
        # p.append(pMean.GetValue(j))
        # arc.append(x.GetValue(j))
        # #print T1Mean.GetTuple(m)[0]
    # #print p
    # try:
        # sT=[x + y for x, y in zip(sT, T)]
        # sP=[x + y for x, y in zip(sP, p)]
    #polyData =  pdi.GetBlock(0).GetPointData()
    #print polyData
    #a=[0.0, 0.0,0.0]
    #print pdi.GetBlock(0).GetPoint(1,a)
    #print a
    #print 
    


#Delete(slice1)
#del slice1

