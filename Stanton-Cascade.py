import numpy as np
import os
import sys
#sys.modules[__name__].__dict__.clear()
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
print(os.getcwd())
#path1=np.genfromtxt("/home/03624/ykanani/Repositories/paraview-scripts/localOutputFolder",dtype='str')
#print(path1)
##############################
############# 1.55/100000/0.711*38000/(15.8*(T1Mean-20))

n= 50
zmin=0.117 #0.117
zmax=0.137 #0.137
##FOR 2D case
#zmin=0.01
#zmax=0.02
Pr = 0.711
nu = 1.55/100000.0
alfa = nu/Pr
gradT = 38000.0
Uexit = 15.6
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

location =input("Please enter the location that you are running this script, 0:Office PC, 1:Marvericks :")
if location==0:
	path1="/home/03624/ykanani/tempparaview/"
elif location==1:
	path1='D:/PostProcess/pvout/'
path2=".plt"


path = path1[0] + caseName + "_St_" +str(zmin) + "to" +str(zmax) + "_n" + str(n)  + path2
fout = open(path,'a')
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
    Curves1.SliceType.Origin = [0.1421842037884744, -0.15909159556031227, z]
    disp = Show(Curves1, renderView1)
    Curves1.UpdatePipeline()
    

    pdi = servermanager.Fetch(Curves1)
    
    #print pdi
    ###print pdi.GetBlock(0)
    #print pdi.GetBlock(0).GetBlock(0)
    #print pdi.GetBlock(0).GetBlock(0).GetBlock(0).GetPointData()
    T1Mean = pdi.GetBlock(0).GetBlock(0).GetBlock(0).GetPointData().GetArray("T1Mean")  #T1 or T1Mean
    T2Mean = pdi.GetBlock(0).GetBlock(0).GetBlock(0).GetPointData().GetArray("T2Mean") #T2 or T2Mean
    pMean = pdi.GetBlock(0).GetBlock(0).GetBlock(0).GetPointData().GetArray("pMean")  #p or pMean
    x = pdi.GetBlock(0).GetBlock(0).GetBlock(0).GetPointData().GetArray("arc_length")
    print "here"
    m=T1Mean.GetNumberOfTuples()
    T=[]
    p=[]
    arc=[]
    #print T1Mean.GetTuple(0)[0]
    for j in range(0,m):
        T.append(T1Mean.GetValue(j))
        p.append(pMean.GetValue(j))
        arc.append(x.GetValue(j))
        #print T1Mean.GetTuple(m)[0]
    #print p
    try:
        sT=[x + y for x, y in zip(sT, T)]
        sP=[x + y for x, y in zip(sP, p)]
    except:
        print "except"
        print i
        sT=T
        sP=p
print sT
print T
avgT=[x/n for x in sT]
avgP=[x/n for x in sP]

dT = [x-20.0 for x in avgT]
#print avgP
Pmax=max([float(i) for i in avgP])
#print "Pmax =" + str(max(avgP))
pnon = [(Pmax-x)/(0.5*Uexit**2) for x in avgP]
PmaxIndex= avgP.index(Pmax)
print PmaxIndex
print Pmax
print arc[234]
xshift=arc[PmaxIndex]
print("Xshift =====================",xshift)
#total length = 1.10812
print arc[-1]
#arcNew=[x-0.63 for x in arc] ##decomposed case
arcNew=[x-xshift for x in arc]   #shift stagnation point to zero
#print arc
print arcNew[1], arcNew[-1]
arcNew=[x if x<0.63 else x-arc[-1] for x in arcNew] #rearrange arc between -.46 and 0.6 ##reconstructed case
print arcNew[1], arcNew[-1]
arcNew=[x if x>-0.47812 else x+arc[-1] for x in arcNew] #rearrange arc between -.46 and 0.6 ##reconstructed case
print arcNew[1], arcNew[-1]
#print dT
St = [alfa * gradT / (Uexit*x) for x in dT]
#print St
#print avgT
#print arc
#print type(data)
data=np.c_[arcNew,St, pnon]

datasorted = data[data[:, 0].argsort()]
#np.savetxt(fout,St,fmt="%10.5f ")
np.savetxt(fout,datasorted,fmt="%15.12f %15.12f %15.12f ")
fout.close()
del sT
del sP
del St
del avgT
del avgP
del arc
del arcNew
del T
del data
