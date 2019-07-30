##NOTES: output parameters need to be updates similar to cascade-lati-avg-ver2.py
#help("modules")

#print(sys.path)

import numpy as np
import os



#sys.modules[__name__].__dict__.clear()
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
##############################

##Flow properties
Uinf = 0.48
delta0 = 0.005641656
##Fluid properties
Pr = 0.711
nu = 9.5616E-7
alfa = nu/Pr
W=0.0498
D=0.006


#################
mm = 500
first=1e-9
maxy=0.025#2*delta0# 3*delta0 harttnet case
l=np.logspace(np.log10(first), np.log10(maxy), num=mm,base=10)
l=np.insert(l,0,0)
#print l
nn = 2
xd0=0.006
Smin=0.0001 #-13 harttnet case
Smax=0.19# -1 hartnett case
ylist=np.array([0.15, 0.175, 0.25])*D   #([0 ,.15 ,.15+0.025 ,.15-0.025])
zlist=np.array([0])*D
#nS=(Smax-Smin)/deltaS+1
#SS=np.linspace(Smin,Smax,nS)3


##############################


xi=Smin
yi=0

##############################
##Needs review for heat transfer
gradT = 38000.0
T0=20


#############
##############################
#############
###############################
## find source
sourceName = raw_input("Please enter source name: ")
#sourceName ="dheh2"
#tname =input("Please enter 0:no heat transfer, 1:T, 2:T1 and T2")
#sourceName = "CL2"
#tname =input("Please enter 0 for T1, 1 for T2")
caseName = sourceName # raw_input("Please enter case name: ")
print sourceName
source = FindSource(sourceName + '.foam')
sourceS = FindSource(sourceName + 'S.foam')
print source
# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

plotOverLine1 = PlotOverLine(Input=source,
    Source='High Resolution Line Source')
#start = input("Please enter start: ")
#end = input("Please enter end: ")
#sp=0.031750
# init the 'High Resolution Line Source' selected for 'Source'
d = 0.1
Uj=input("Please enter Uj")


plotOverLine1.Source.Point1 = [Smin, yi+d, 0.0]
plotOverLine1.Source.Point2 = [Smax, yi+d, 0.0]

# Properties modified on plotOverLine1
plotOverLine1.Tolerance = 2.22044604925031e-16
# get active view
renderView1 = GetActiveViewOrCreate('RenderView')



#####################################surface
slice2 = Slice(Input=sourceS)
slice2.SliceType = 'Plane'
slice2.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice2.SliceType.Origin = [xi, yi, 0]
slice2.SliceType.Normal = [0.0,0.0, 1]


############################################
##############################
##############################
##############################

####### initializing files #######
#path1="/home/03624/ykanani/tempparaview/"
path1="E:/PostProcess/pvout/"
path2=".plt"


path = path1 + caseName + "_Xprofiles_NoAvg"  + path2
foutv2 = open(path,'a')



try:
        del y
except NameError:
        pass

try:
        del UprofAvg
except NameError:
        pass

try:
        del TprofAvg
except NameError:
        pass
index=1

for dz in zlist:
    for dy in ylist:
        print "Current Location Y/D = " + str(dy/D)
        print "Current Location Z/D = " + str(dz/D)
        yi =  dy
        zi= dz
        foutv2.write( "VARIABLES= X/D Y/D Z/D T1 T2 T1rms T2rms Ux/Uj Uy/Uj Uz/Uj UU VV WW UV VW UW" + " \n" )
        foutv2.write( "ZONE T=\" " + sourceName + "_Y=" + str(dy/D) + "_Z=" + str(dz/D) + "\" \n" )
        T1Avg=[]
        T2Avg=[]
        T1rmsAvg=[]
        T2rmsAvg=[]
            
        UAvg=[]
        VAvg=[]
        WAvg=[]
            
        UrmsAvg=[]
        VrmsAvg=[]
        WrmsAvg=[]
        UVAvg=[]
        VWAvg=[]
        UWAvg=[]
        x=[]

            
        plotOverLine1.Source.Point1 = [Smin, yi, zi]
        plotOverLine1.Source.Point2 = [Smax, yi, zi]

        slice2.SliceType.Origin = [0, 0, zi]
        
        if dy ==0 :
            pdi = servermanager.Fetch(slice2)
            print(pdi.GetBlock(0).GetBlock(0))
            polydata=pdi.GetBlock(0).GetBlock(0)
        else:
            pdi = servermanager.Fetch(plotOverLine1)
            polydata=pdi

            
        T1Mean = polydata.GetPointData().GetArray("T1Mean")  #T1 or T1Mean
        T2Mean = polydata.GetPointData().GetArray("T2Mean") #T2 or T2Mean
        T1Prime2Mean = polydata.GetPointData().GetArray("T1Prime2Mean")  #T1Prime2Mean
        T2Prime2Mean = polydata.GetPointData().GetArray("T2Prime2Mean")  #T2Prime2Mean                   
        UMean = polydata.GetPointData().GetArray("UMean")  #p or pMean
        UPrime2Mean = polydata.GetPointData().GetArray("UPrime2Mean")  #p or pMean , UPrime2Mean
        arc = polydata.GetPointData().GetArray("arc_length")
        m=T1Mean.GetNumberOfTuples()
        T1=[]
        T2=[]
        T1ms=[]
        T2ms=[]
        U=[]
        V=[]
        W=[]
        UU=[]
        VV=[]
        WW=[]
        UV=[]
        VW=[]
        UW=[]
        X=[]
        #print T1Mean.GetTuple(0)[0]
        for j in range(0,m):
            T1.append(T1Mean.GetValue(j))
            T2.append(T2Mean.GetValue(j))
            T1ms.append(T1Prime2Mean.GetValue(j))
            T2ms.append(T2Prime2Mean.GetValue(j))
            
            U.append(UMean.GetTuple(j)[0]/Uj)
            V.append(UMean.GetTuple(j)[1]/Uj)
            W.append(UMean.GetTuple(j)[2]/Uj)
            print(T1Mean)
            print(UPrime2Mean.GetTuple(j))
            UU.append(UPrime2Mean.GetTuple(j)[0])
            VV.append(UPrime2Mean.GetTuple(j)[1])
            WW.append(UPrime2Mean.GetTuple(j)[2])
            UV.append(UPrime2Mean.GetTuple(j)[3])
            VW.append(UPrime2Mean.GetTuple(j)[4])
            UW.append(UPrime2Mean.GetTuple(j)[5])
            X.append(arc.GetValue(j)/D)

                
        dyArray = np.ones(len(X))*dy/D
        dzArray = np.ones(len(X))*dz/D
        np.savetxt(foutv2,np.c_[X,dyArray,dzArray,T1,T2,T1ms,T2ms,U,V,W,UU,VV,WW,UV,VW,UW],\
                        fmt="%15.12f "*16)                                                                                                                                                                                                                                                            

foutv2.close()
del foutv2


Delete(plotOverLine1)
del plotOverLine1



Delete(slice2)
del slice2

