import numpy as np
import os
#import sys
#sys.modules[__name__].__dict__.clear()
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

##############################
############# 1.55/100000/0.711*38000/(15.8*(T1Mean-20))

n= 300
zmin=0.1
zmax=4.9
##FOR 2D case
#zmin=0.01
#zmax=0.02
Pr = 0.711
nu = 0.00013895
alfa = nu/Pr
gradT=1
Uinf=1

#############
###############################
## find source
sourceName = raw_input("Please enter source name: ")
#sourceName = "CL2"
#tname =input("Please enter 0 for T1, 1 for T2")
caseName = sourceName # raw_input("Please enter case name: ")
print sourceName
source = FindSource(sourceName + '.foam')
print source
# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# create a new 'Plot On Intersection Curves'
Curves1 = PlotOnIntersectionCurves(Input=source)
Curves1.SliceType = 'Plane'

# init the 'Plane' selected for 'SliceType'
Curves1.SliceType.Origin = [0.1421842037884744, -0.15909159556031227, 0.12700000405311584]
Curves1.SliceType.Normal = [0, 0, 1]


disp = Show(Curves1, renderView1)


#path1="/home/03624/ykanani/tempparaview/"
path1='D:/PostProcess/pvout/'
path2=".plt"


path = path1 + caseName + "_eta_st"  + path2

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).

for i in range(0,n):
    print i
    if (n==1):
        z=zmin
    else:
        z = (zmax-zmin)/(n-1)*i+zmin
    
    
    print z
    Curves1.SliceType.Origin = [0, 0, z]
    disp = Show(Curves1, renderView1)
    Curves1.UpdatePipeline()
    

    pdi = servermanager.Fetch(Curves1)
    
    #print pdi
    #print pdi.GetBlock(0)
    #print pdi.GetBlock(0).GetBlock(0)
    #print pdi.GetBlock(0).GetBlock(0).GetBlock(0).GetPointData()
    T1Mean = pdi.GetBlock(0).GetBlock(0).GetBlock(0).GetPointData().GetArray("T1Mean")  #T1 or T1Mean
    T2Mean = pdi.GetBlock(0).GetBlock(0).GetBlock(0).GetPointData().GetArray("T1Mean") #T2 or T2Mean
    pMean = pdi.GetBlock(0).GetBlock(0).GetBlock(0).GetPointData().GetArray("T1Mean")  #p or pMean
    xarc = pdi.GetBlock(0).GetBlock(0).GetBlock(0).GetPointData().GetArray("arc_length")

    m=T1Mean.GetNumberOfTuples()
    T1=[]
    T2=[]
    p=[]
    X=[]
    #print T1Mean.GetTuple(0)[0]
    for j in range(0,m):
        T1.append(T1Mean.GetValue(j))
        T2.append(T2Mean.GetValue(j))
        p.append(pMean.GetValue(j))
        X.append(xarc.GetValue(j))
        #print T1Mean.GetTuple(m)[0]
    #print p
    try:
        sT1=[x + y for x, y in zip(sT1, T1)]
        sT2=[x + y for x, y in zip(sT2, T2)]
        sP=[x + y for x, y in zip(sP, p)]
    except:
        print "except"
        print i
        sT1=T1
        sT2=T2
        sP=p
print sT1
print T1
avgT1=[x/n for x in sT1]
avgT2=[x/n for x in sT2]
avgP=[x/n for x in sP]


Pmax=max(avgP)
#print "Pmax =" + str(max(avgP))
pnon = [(Pmax-x)/(0.5*Uinf**2) for x in avgP]
PmaxIndex= avgP.index(Pmax)
print PmaxIndex



print avgT1
print avgT2


dT2 = [x-y for x, y in zip(avgT2, avgT1)]
St = [alfa * gradT / (Uinf*x) for x in dT2]
eta = avgT1

xp=[x+4.065040634154 for x in X] 
xms=[(x)/0.28 for x in X] 
xpms=[(x+4.065040634154)/0.28 for x in X] 
#print St
#print avgT
#print arc
#print type(data)
data=np.c_[xp,X,xpms,xms,eta,St, pnon]

datasorted = data[data[:, 0].argsort()]

fout = open(path,'a')
fout.write( "VARIABLES= xp x xpms xms eta st pnon" + " \n" )
fout.write( "ZONE T=\" " + sourceName + "_plate" + "\" \n" )
np.savetxt(fout,datasorted,fmt="%15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f")
fout.close()
del sT1
del sT2
del sP
del St
del eta
del avgT1
del avgT2
del avgP
del x
del T1
del T2
del data
