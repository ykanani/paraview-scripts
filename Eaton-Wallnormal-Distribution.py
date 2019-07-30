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
#Uinf = 0.48
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
Smin=-10*D #-13 harttnet case
Smax=0.2# -1 hartnett case
deltaS=5*D
#SS=np.arange(deltaS,Smax,deltaS)
SS=np.array([-8,-7,-6,-5,-4,-3,-2,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,25,30])*D

zmin=-0.5*D
zmax=+0.6*D
zlist=np.arange(zmin,zmax,0.5*D)

ymax=4*D

#nS=(Smax-Smin)/deltaS+1
#SS=np.linspace(Smin,Smax,nS)
print SS

##############################

n= 5
#delta=0.02

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
#sourceName ="dheh2"
tname =input("Please enter 0:no heat transfer, 1:T, 2:T1 and T2")
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
Uj=input("Please enter Uj")


plotOverLine1.Source.Point1 = [xi, 0, zmin]
plotOverLine1.Source.Point2 = [xi, ymax, zmin]

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')





############################################
##############################
##############################
##############################

####### initializing files #######
#path1="/home/03624/ykanani/tempparaview/"
path1="E:/PostProcess/pvout/"
path2=".plt"


###saving PROFILES
path = path1 + caseName + "_YProfiles" + path2
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
    for dz in zlist:
            print('X/D = ',str(dx/D),'Z/D = ',str(dz/D))
            foutv2.write( "VARIABLES= X/D Y/D Z/D T1 T2 T1rms T2rms Ux/Uj Uy/Uj Uz/Uj Urms Vrms Wrms UV VW UW" + " \n" )
            foutv2.write( "ZONE T=\" " + "case1" + "_X=" + str(dx/D) + "_Z=" + str(dz/D) + "\" \n" )
            #print j
            #print "d = " + str(d)
            plotOverLine1.Source.Point1 = [xi, 0, dz]
            plotOverLine1.Source.Point2 = [xi, ymax, dz]

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
            y=[]
            
            #print T1Mean.GetTuple(0)[0]
            #print(UPrime2Mean.GetTuple(0))
            for j in range(0,m):
                T1.append(T1Mean.GetValue(j))
                T2.append(T2Mean.GetValue(j))
                T1rms.append(np.sqrt(T1Prime2Mean.GetValue(j)))
                T2rms.append(np.sqrt(T2Prime2Mean.GetValue(j)))
                
                U.append(UMean.GetTuple(j)[0]/Uj)
                V.append(UMean.GetTuple(j)[1]/Uj)
                W.append(UMean.GetTuple(j)[2]/Uj)
                
                Urms.append(np.sqrt(UPrime2Mean.GetTuple(j)[0]))
                Vrms.append(np.sqrt(UPrime2Mean.GetTuple(j)[1]))
                Wrms.append(np.sqrt(UPrime2Mean.GetTuple(j)[2]))
                UV.append(UPrime2Mean.GetTuple(j)[3])
                VW.append(UPrime2Mean.GetTuple(j)[4])
                UW.append(UPrime2Mean.GetTuple(j)[5])
                
                y.append((arc.GetValue(j))/D)
                
            
            #print p
            print len(y)
            print len(T1)
            dxArray = np.ones(len(y))*dx/D
            dzArray = np.ones(len(y))*dz/D
            np.savetxt(foutv2,np.c_[dxArray,y,dzArray,T1,T2,T1rms,T2rms,U,V,W,Urms,Vrms,Wrms,UV,VW,UW],\
                    fmt="%15.12f "*16)

    

del SS
foutv2.close()
del foutv2


Delete(plotOverLine1)
del plotOverLine1

                                                                                                                                                                                                                                                                                                        
