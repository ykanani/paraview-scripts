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
Smin=0.0001 #-13 harttnet case
Smax=0.32# -1 hartnett case
deltaS=0.01
SS=np.arange(deltaS,Smax,deltaS)
SS=np.insert(SS,0,Smin)
ymad=2*D
ylist=np.array([0,0.15,0.25,0.35])
#nS=(Smax-Smin)/deltaS+1
#SS=np.linspace(Smin,Smax,nS)3
print SS

##############################

n= 5
#delta=0.02
zmin=-2*D
zmax=+2*D
sp=zmax-zmin  #spanwise distance

print sp
#############################

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
tname =input("Please enter 0:no heat transfer, 1:T, 2:T1 and T2")
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



plotOverLine1.Source.Point1 = [xi, yi+d, zmin]
plotOverLine1.Source.Point2 = [xi, yi+d, zmax]

# Properties modified on plotOverLine1
plotOverLine1.Tolerance = 2.22044604925031e-16
# get active view
renderView1 = GetActiveViewOrCreate('RenderView')



# create a new 'Integrate Variables'
integrateVariables1 = IntegrateVariables(Input=plotOverLine1)

#####################################surface
clip1 = Clip(Input=sourceS)
clip1.ClipType = 'Plane'
clip1.ClipType.Origin = [0, 0, zmin]
clip1.ClipType.Normal = [0,0,1]
clip1.Invert=0
clip2 = Clip(Input=clip1)
clip2.ClipType = 'Plane'
clip2.ClipType.Origin = [0, 0, zmax]
clip2.ClipType.Normal = [0,0,-1]
clip2.Invert=0
# create a new 'Slice' on the surface
slice2 = Slice(Input=clip2)
slice2.SliceType = 'Plane'
slice2.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice2.SliceType.Origin = [xi, yi, 0]
slice2.SliceType.Normal = [1,0.0, 0.0]

integrateVariables2 = IntegrateVariables(Input=slice2)

############################################
##############################
##############################
##############################

####### initializing files #######
#path1="/home/03624/ykanani/tempparaview/"
path1="E:/PostProcess/pvout/"
path2=".plt"


path = path1 + caseName + "_Xprofiles"  + path2
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

for dy in ylist:
    print "Current Location Y/D + " + str(y)
    yi =  dy
    foutv2.write( "VARIABLES= X/D Y/D T1 T2 T1rms T2rms Ux/Uj Uy/Uj Uz/Uj Urms Vrms Wrms UV VW UW" + " \n" )
    foutv2.write( "ZONE T=\" " + "case1" + "_Y=" + str(dy/D) + "\" \n" )
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
    for dS in SS:
            
        xi=dS
        
        plotOverLine1.Source.Point1 = [xi, yi, zmin]
        plotOverLine1.Source.Point2 = [xi, yi, zmax]

        slice2.SliceType.Origin = [xi, yi, 0]

        if dy ==0 :
            pdi = servermanager.Fetch(integrateVariables2)
            polydata=pdi
        else:
            pdi = servermanager.Fetch(integrateVariables1)
            polydata=pdi

        
        T1Mean = polydata.GetPointData().GetArray("T1Mean")  #T1 or T1Mean
        T2Mean = polydata.GetPointData().GetArray("T2Mean") #T2 or T2Mean
        T1Prime2Mean = polydata.GetPointData().GetArray("T1Prime2Mean")  
        T2Prime2Mean = polydata.GetPointData().GetArray("T2Prime2Mean")                     
        UMean = polydata.GetPointData().GetArray("UMean")  #p or pMean
        UPrime2Mean = polydata.GetPointData().GetArray("UPrime2Mean")  #p or pMean
        
        T1Avg.append(T1Mean.GetValue(0))
        T2Avg.append(T2Mean.GetValue(0))
        T1rmsAvg.append(np.sqrt(T1Prime2Mean.GetValue(0))/sp)
        T2rmsAvg.append(np.sqrt(T2Prime2Mean.GetValue(0))/sp)
            
        UAvg.append(UMean.GetTuple(0)[0]/sp)
        VAvg.append(UMean.GetTuple(0)[1]/sp)
        WAvg.append(UMean.GetTuple(0)[2]/sp)
            
        UrmsAvg.append(np.sqrt(UPrime2Mean.GetTuple(0)[0])/sp)
        VrmsAvg.append(np.sqrt(UPrime2Mean.GetTuple(0)[1])/sp)
        WrmsAvg.append(np.sqrt(UPrime2Mean.GetTuple(0)[2])/sp)
        UVAvg.append(UPrime2Mean.GetTuple(0)[3]/sp)
        VWAvg.append(UPrime2Mean.GetTuple(0)[4]/sp)
        UWAvg.append(UPrime2Mean.GetTuple(0)[5]/sp)
        
        
        x.append(xi/D)
            
    dyArray = np.ones(len(x))*dy/D
    np.savetxt(foutv2,np.c_[dxArray,dyArray,T1,T2,T1rms,T2rms,U,V,W,Urms,Vrms,Wrms,UV,VW,UW],\
                    fmt="%15.12f "*15)                                                                                                                                                                                                                                                            
