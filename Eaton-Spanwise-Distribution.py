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

#print l
nn = 2
Smin=0 #-13 harttnet case
Smax=0.2# -1 hartnett case
deltaS=5*D
SS=np.arange(deltaS,Smax,deltaS)
#ylist=np.arange(0.0,0.45*D,0.15*D)
ylist=np.array([0.25*D])
#nS=(Smax-Smin)/deltaS+1
#SS=np.linspace(Smin,Smax,nS)
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
#sourceName = raw_input("Please enter source name: ")
sourceName ="dheh2"
#tname =input("Please enter 0:no heat transfer, 1:T, 2:T1 and T2")
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
d = 0.01



plotOverLine1.Source.Point1 = [xi, yi+d, zmin]
plotOverLine1.Source.Point2 = [xi, yi+d, zmax]

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')



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
slice2 = PlotOnIntersectionCurves(Input=clip2)
slice2.SliceType = 'Plane'
# init the 'Plane' selected for 'SliceType'
slice2.SliceType.Origin = [xi, yi, 0]
slice2.SliceType.Normal = [1,0.0, 0.0]

############################################
##############################
##############################
##############################

####### initializing files #######
#path1="/home/03624/ykanani/tempparaview/"
path1="E:/PostProcess/pvout/"
path2=".plt"


###saving PROFILES
path = path1 + caseName + "_ZProfiles" + path2
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

for dx in SS:


    
    xi = dx
    yi = 0

    
    

    
    #foutv2.write( "ZONE T=\" " + sourceName + "_" + str(dS) + "\" \n" )
    j=0
    for dy in ylist:
            print('X/D = ',str(dx/D),'Y/D = ',str(dy/D))
            foutv2.write( "VARIABLES= x y z T1 T2 T1rms T2rms U V W Urms Vrms Wrms UV VW UW" + " \n" )
            foutv2.write( "ZONE T=\" " + "case1" + "_X=" + str(dx/D) + "_Y=" + str(dy/D) + "\" \n" )
            #print j
            #print "d = " + str(d)
            plotOverLine1.Source.Point1 = [xi, yi+dy, zmin]
            plotOverLine1.Source.Point2 = [xi, yi+dy, zmax]
            slice2.SliceType.Origin = [xi, 0, 0]

            if dy ==0 :
                pdi = servermanager.Fetch(slice2)
                polydata=pdi.GetBlock(0).GetBlock(0).GetBlock(0)
            else:
                pdi = servermanager.Fetch(plotOverLine1)
                polydata=pdi
                
                
            #print pdi
            
            T1Mean = polydata.GetPointData().GetArray("T1Mean")  #T1 or T1Mean
            T2Mean = polydata.GetPointData().GetArray("T2Mean") #T2 or T2Mean
            T1Prime2Mean = polydata.GetPointData().GetArray("T1Prime2Mean")  
            T2Prime2Mean = polydata.GetPointData().GetArray("T2Prime2Mean")                     
            UMean = polydata.GetPointData().GetArray("UMean")  #p or pMean
            UPrime2Mean = polydata.GetPointData().GetArray("UPrime2Mean")  #p or pMean
            
            arc = polydata.GetPointData().GetArray("arc_length")
            m=T1Mean.GetNumberOfTuples()
            T1=[]
            T2=[]
            T1rms=[]
            T2rms=[]
            U=[]
            V=[]
            W=[]
            Urms=[]
            Vrms=[]
            Wrms=[]
            UV=[]
            VW=[]
            UW=[]
            z=[]
            
            #print T1Mean.GetTuple(0)[0]
            #print(UPrime2Mean.GetTuple(0))
            for j in range(0,m):
                T1.append(T1Mean.GetValue(j))
                T2.append(T2Mean.GetValue(j))
                T1rms.append(np.sqrt(T1Prime2Mean.GetValue(j)))
                T2rms.append(np.sqrt(T2Prime2Mean.GetValue(j)))
                
                U.append(UMean.GetTuple(j)[0])
                V.append(UMean.GetTuple(j)[1])
                W.append(UMean.GetTuple(j)[2])
                
                Urms.append(np.sqrt(UPrime2Mean.GetTuple(j)[0]))
                Vrms.append(np.sqrt(UPrime2Mean.GetTuple(j)[1]))
                Wrms.append(np.sqrt(UPrime2Mean.GetTuple(j)[2]))
                UV.append(UPrime2Mean.GetTuple(j)[3])
                VW.append(UPrime2Mean.GetTuple(j)[4])
                UW.append(UPrime2Mean.GetTuple(j)[5])
                
                z.append((arc.GetValue(j)-zmax)/D)
                #print T1Mean.GetTuple(m)[0]
            
            #print p
            dxArray = np.ones(len(z))*dx/D
            dyArray = np.ones(len(z))*dy/D
            np.savetxt(foutv2,np.c_[dxArray,dyArray,z,T1,T2,T1rms,T2rms,U,V,W,Urms,Vrms,Wrms,UV,VW,UW],\
                    fmt="%15.12f "*16)

    

del SS
foutv2.close()
del foutv2

Delete(slice2)
del slice2
Delete(clip2)
del clip2
Delete(clip1)
del clip1

Delete(plotOverLine1)
del plotOverLine1

                                                                                                                                                                                                                                                                                                        
