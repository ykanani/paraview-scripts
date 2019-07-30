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
Uinlet = 0.52
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
deltaS=0.006
SS=np.arange(deltaS,Smax,deltaS)
SS=np.insert(SS,0,Smin)

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

# create a new 'Slice'
slice1 = Slice(Input=source)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [xi, yi, 0]
slice1.SliceType.Normal = [1,0, 0.0]




# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

plotOverLine1 = PlotOverLine(Input=slice1,
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
path1='E:/PostProcess/pvout/'
path2=".plt"


path = path1 + caseName + "_BL_All"  + path2
boundarylayer = open(path,'a')
boundarylayer.write( "VARIABLES= dS Ue pS  St1 St2 T1S T2S T1E T2E delta delta95 delta99 delta999 \
                     deltas theta H G \
                     deltaT1 deltaT195 deltaT199 deltaT1999 \
                     deltaT2 deltaT295 deltaT299 deltaT2999 \
                     deltaT1W deltaT2W \
                     cf1 cf2 utau R UUEdge VVEdge WWEdge" + " \n")
boundarylayer.write( "ZONE T=\" " + sourceName + "_BL_ALL" + "\" \n" )

path = path1 + caseName + "_MAX_All"  + path2
Maxvalues = open(path,'a')
Maxvalues.write( "VARIABLES= dS Ue  T1S T2S T1E T2E UUMax VVMax WWMax UVMax UWMax VWMax UVMaxNon,\
UUEdge, VVEdge,WWEdge, UVEdge, UWEdge, VWEdge, \
TT1Max TT2Max yUUMax yVVMax yWWMax yUVMax yUWMax yVWMax yTT1Max yTT2Max delta delta95 delta99 delta999" + " \n" )
Maxvalues.write( "ZONE T=\" " + sourceName + "_MAX_ALL" + "\" \n" )

###saving PROFILES
path = path1 + caseName + "_Profiles" + path2
Profiles = open(path,'a')

ri = 0 # dummy

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

def intsum(x,dx):
        #print "in intsum ..........."
        #print x
        #print dx
        lenx = len(x)
        #print lenx
        x0 = x[0:lenx-1]
        x1 = x[1:lenx]
        xAve = [sum(x)/2.0 for x in zip(x0,x1)]
        dx0 = dx[0:lenx-1]
        dx1 = dx[1:lenx]
        dxAve = [xx1 - xx0 for xx1 , xx0 in zip(dx1,dx0)]
        
        return np.inner(xAve,dxAve)

for dS in SS:

    # restarting main profile variables
    y=[]
    pprofAvg=[]
    
    UprofAvg=[]
    VprofAvg=[]
    WprofAvg=[]
    
    UVprofAvg=[]
    VWprofAvg=[]
    UWprofAvg=[]
    
    UUprofAvg=[]
    VVprofAvg=[]
    WWprofAvg=[]
    
    gradUprgAvg=[]
    
    nuSgsprofAvg=[]
    
    T1profAvg=[]
    T2profAvg=[]
    
    TT1profAvg=[]
    TT2profAvg=[]
    
    VT1profAvg=[]
    VT2profAvg=[]
    
    gradT1profAvg=[]
    gradT2profAvg=[]
    Sa =  dS
    print "Current Location + " + str(Sa)
    xi = Sa
    yi = 0

    # modify slice1
    slice1.SliceType.Origin = [xi, yi, 0]
    slice1.SliceType.Normal = [1,0.0, 0.0]
    # modify slice 2
    slice2.SliceType.Origin = [xi, yi, 0]
    slice2.SliceType.Normal = [1,0, 0.0]

    Profiles.write( "VARIABLES= dS Ue utau delta delta95 delta99 delta999 deltas theta \
                                T1E T2E T1S T2S \
                                pS \
                                deltaT1W deltaT2W \
                                delta95T1 delta95T2 \
                                delta99T1 delta99T2 \
                                y \
                                p \
                                U V W \
                                UU VV WW UV UW VW \
                                gradU \
                                nuSgs \
                                T1 T2 \
                                ThetaT1 ThetaT2 \
                                TT1 TT2 \
                                VT1 VT2 \
                                UT1 UT2 \
                                Tu \
                                gradT1 gradT2 \
                                " + " \n" )
    Profiles.write( "ZONE T=\" " + sourceName + "_" + str(index) + "_" +  str(round(dS/D,2)) + "\" \n" )
    j=0
    
    for d in l:
        
        #d = maxy/(mm-1)*j
        #d=1e-6*(1-1.01**j)/(1-1.01)
        #print j
        #print "d = " +  str(d)
        plotOverLine1.Source.Point1 = [xi, yi+d, zmin]
        plotOverLine1.Source.Point2 = [xi, yi+d, zmax]



        pdi = servermanager.Fetch(integrateVariables1)
        pdiS = servermanager.Fetch(integrateVariables2)
        
        #print pdi
        if j ==0 :
            try:
                wgum = pdiS.GetPointData().GetArray("wallGradUMean")
                pMeanS = pdiS.GetPointData().GetArray("pMean")
                if tname==0:
                    T1S=0
                    T2S=0
                    TT1S=0
                    TT2S=0
                    gradT2S=0
                elif tname==1:
                    T1MeanS = pdiS.GetPointData().GetArray("TMean")
                    T2MeanS = pdiS.GetPointData().GetArray("TMean")
                    T1Prime2MeanS = pdiS.GetPointData().GetArray("TPrime2Mean")
                    T2Prime2MeanS = pdiS.GetPointData().GetArray("TPrime2Mean")
                    gradT1MeanS = pdiS.GetPointData().GetArray("gradTMeantn")
                    gradT2MeanS = pdiS.GetPointData().GetArray("gradTMeantn")
                elif tname==1:
                    T1MeanS = pdiS.GetPointData().GetArray("T1Mean")
                    T2MeanS = pdiS.GetPointData().GetArray("T1Mean")
                    T1Prime2MeanS = pdiS.GetPointData().GetArray("T1Prime2Mean")
                    T2Prime2MeanS = pdiS.GetPointData().GetArray("T1Prime2Mean")
                    gradT1MeanS = pdiS.GetPointData().GetArray("gradT1Meantn")
                    gradT2MeanS = pdiS.GetPointData().GetArray("gradT1Meantn")
                elif tname==2:
                    T1MeanS = pdiS.GetPointData().GetArray("T1Mean")
                    T2MeanS = pdiS.GetPointData().GetArray("T2Mean")
                    T1Prime2MeanS = pdiS.GetPointData().GetArray("T1Prime2Mean")
                    T2Prime2MeanS = pdiS.GetPointData().GetArray("T2Prime2Mean")
                    gradT1MeanS = pdiS.GetPointData().GetArray("gradT1Meantn")
                    gradT2MeanS = pdiS.GetPointData().GetArray("gradT2Meantn")
                    
                wgu = np.sqrt(wgum.GetTuple(0)[0]**2+wgum.GetTuple(0)[1]**2)/sp
                pS = pMeanS.GetTuple(0)[0]/sp
                if tname!=0:
                    T1S = T1MeanS.GetTuple(0)[0]/sp
                    T2S = T2MeanS.GetTuple(0)[0]/sp
                    TT1S = T1Prime2MeanS.GetTuple(0)[0]/sp
                    TT2S = T1Prime2MeanS.GetTuple(0)[0]/sp
                    gradT1S = gradT1MeanS.GetTuple(0)[0]/sp
                    gradT2S = gradT2MeanS.GetTuple(0)[0]/sp
                    print "GradT1S=" + str(gradT1S)
                    print "GradT2S=" + str(gradT2S)
                
                print "WGU =" + str(wgu)
            except:
                #not on the wall
                print "No wall or WallgradU is not available!!! skipping this dS/D =" + str(dS/D)
                wgu = -1
                pS = 0
                if tname==1:
                    T1MeanS = pdi.GetPointData().GetArray("TMean")
                    T2MeanS = pdi.GetPointData().GetArray("TMean")
                    T1Prime2MeanS = pdi.GetPointData().GetArray("TPrime2Mean")
                    T2Prime2MeanS = pdi.GetPointData().GetArray("TPrime2Mean")
                    gradT1MeanS = pdi.GetPointData().GetArray("gradTMean")
                    gradT2MeanS = pdi.GetPointData().GetArray("gradTMean")
                elif tname==2:
                    T1MeanS = pdi.GetPointData().GetArray("T1Mean")
                    T2MeanS = pdi.GetPointData().GetArray("T2Mean")
                    T1Prime2MeanS = pdi.GetPointData().GetArray("T1Prime2Mean")
                    T2Prime2MeanS = pdi.GetPointData().GetArray("T2Prime2Mean")
                    gradT1MeanS = pdi.GetPointData().GetArray("gradT1Mean")
                    gradT2MeanS = pdi.GetPointData().GetArray("gradT1Mean")
                if tname!=0:
                    T1S = T1MeanS.GetTuple(0)[0]/sp
                    T2S = T2MeanS.GetTuple(0)[0]/sp
                    TT1S = T1Prime2MeanS.GetTuple(0)[0]/sp
                    TT2S = T1Prime2MeanS.GetTuple(0)[0]/sp
        
        #reading variables
        pMean = pdi.GetPointData().GetArray("pMean")
        UMean = pdi.GetPointData().GetArray("UMean")
        Ums = pdi.GetPointData().GetArray("UPrime2Mean")
        gradUMean = pdi.GetPointData().GetArray("gradUMean")   
        nuSgs = pdi.GetPointData().GetArray("nut")   ##############################FOR LAMINAR I CHANGED NUT TO PMEAN !!!!
        if tname==1:
            T1Mean = pdi.GetPointData().GetArray("TMean")
            T2Mean = pdi.GetPointData().GetArray("TMean")
            
            T1Prime2Mean = pdi.GetPointData().GetArray("TPrime2Mean")
            T2Prime2Mean = pdi.GetPointData().GetArray("TPrime2Mean")
            
            UT1Mean = pdi.GetPointData().GetArray("utAvg")
            UT2Mean = pdi.GetPointData().GetArray("utAvg")
            
            gradT1Mean = pdi.GetPointData().GetArray("gradTMean")
            gradT2Mean = pdi.GetPointData().GetArray("gradTMean")
        elif tname==2:
            T1Mean = pdi.GetPointData().GetArray("T1Mean")
            T2Mean = pdi.GetPointData().GetArray("T1Mean")
            
            T1Prime2Mean = pdi.GetPointData().GetArray("T1Prime2Mean")
            T2Prime2Mean = pdi.GetPointData().GetArray("T1Prime2Mean")
            
            UT1Mean = pdi.GetPointData().GetArray("ut1Avg")
            UT2Mean = pdi.GetPointData().GetArray("ut1Avg")
            
            gradT1Mean = pdi.GetPointData().GetArray("gradT1Mean")
            gradT2Mean = pdi.GetPointData().GetArray("gradT1Mean")
        elif tname==3:
            T1Mean = pdi.GetPointData().GetArray("T1Mean")
            T2Mean = pdi.GetPointData().GetArray("T2Mean")
            
            T1Prime2Mean = pdi.GetPointData().GetArray("T1Prime2Mean")
            T2Prime2Mean = pdi.GetPointData().GetArray("T2Prime2Mean")
            
            UT1Mean = pdi.GetPointData().GetArray("ut1Avg")
            UT2Mean = pdi.GetPointData().GetArray("ut2Avg")
            
            gradT1Mean = pdi.GetPointData().GetArray("gradT1Mean")
            gradT2Mean = pdi.GetPointData().GetArray("gradT2Mean")    
        if d==0:
            #print "if d =0"
            y = [0.0]
            pprofAvg=[pS]
            
            UprofAvg= [0.0]
            VprofAvg=[0.0]
            WprofAvg=[0.0]
            
            UVprofAvg=[0.0]
            UWprofAvg=[0.0]
            VWprofAvg=[0.0]
            
            UUprofAvg=[0.0]
            VVprofAvg=[0.0]
            WWprofAvg=[0.0]
            
            gradUprofAvg=[wgu]
            
            nuSgsprofAvg=[0.0]
            
            T1profAvg= [T1S]
            T2profAvg= [T2S]
            
            TT1profAvg=[TT1S]
            TT2profAvg=[TT2S]
            
            VT1profAvg=[0.0]
            VT2profAvg=[0.0]
            
            UT1profAvg=[0.0]
            UT2profAvg=[0.0]
            
            TuprofAvg=[0.0]
            
            gradT1profAvg=[-38000]
            gradT2profAvg=[gradT2S]
        elif  np.isfinite(UMean.GetTuple(0)[0]/sp):
            
            #print "d="+str(d)
            y.append(d)
            pprofAvg.append(pMean.GetTuple(0)[0]/sp)
            
            UprofAvg.append(UMean.GetTuple(0)[0]/sp)
            VprofAvg.append(UMean.GetTuple(0)[1]/sp)
            WprofAvg.append(UMean.GetTuple(0)[2]/sp)
            
            UVprofAvg.append(Ums.GetTuple(0)[3]/sp)
            VWprofAvg.append(Ums.GetTuple(0)[4]/sp)
            UWprofAvg.append(Ums.GetTuple(0)[5]/sp)
            
            UUprofAvg.append(Ums.GetTuple(0)[0]/sp)
            VVprofAvg.append(Ums.GetTuple(0)[1]/sp)
            WWprofAvg.append(Ums.GetTuple(0)[2]/sp)
            
            gradUprofAvg.append(0)#gradUMean.GetTuple(0)[3]/sp)
            
            nuSgsprofAvg.append(nuSgs.GetTuple(0)[0]/sp)
            if tname!=0:
                T1profAvg.append(T1Mean.GetTuple(0)[0]/sp)
                T2profAvg.append(T2Mean.GetTuple(0)[0]/sp)
                
                TT1profAvg.append(T1Prime2Mean.GetTuple(0)[0]/sp)
                TT2profAvg.append(T2Prime2Mean.GetTuple(0)[0]/sp)
                
                VT1profAvg.append(UT1Mean.GetTuple(0)[1]/sp)
                VT2profAvg.append(UT2Mean.GetTuple(0)[1]/sp)
                
                UT1profAvg.append(UT1Mean.GetTuple(0)[0]/sp)
                UT2profAvg.append(UT2Mean.GetTuple(0)[0]/sp)
                
                TuprofAvg.append(np.sqrt((Ums.GetTuple(0)[0]/sp+Ums.GetTuple(0)[1]/sp+Ums.GetTuple(0)[2]/sp)/3))
                
                gradT1profAvg.append(gradT1Mean.GetTuple(0)[1]/sp)
                gradT2profAvg.append(gradT2Mean.GetTuple(0)[1]/sp)
            else:
                T1profAvg.append(0)
                T2profAvg.append(0)
                
                TT1profAvg.append(0)
                TT2profAvg.append(0)
                
                VT1profAvg.append(0)
                VT2profAvg.append(0)
                
                UT1profAvg.append(0)
                UT2profAvg.append(0)
                
                TuprofAvg.append(0)
                
                gradT1profAvg.append(0)
                gradT2profAvg.append(0)
        j=j+1
    print UprofAvg
    print all([UprofAvg[i] <=0 for i in range(len(UprofAvg)) ])
    UprofAvgAbs=[]
    UprofAvgAbs = map(abs,UprofAvg)
    Ue = max(UprofAvgAbs)
    UeIndex = UprofAvgAbs.index(Ue)
    U95Index = [i for i in range(len(UprofAvg)) if UprofAvgAbs[i] >= 0.95*Ue][0]
    U99Index = [i for i in range(len(UprofAvg)) if UprofAvgAbs[i] >= 0.99*Ue][0]
    U999Index = [i for i in range(len(UprofAvg)) if UprofAvgAbs[i] >= 0.999*Ue][0]
    
    U95= [x for x in UprofAvgAbs if x >= 0.95*Ue][0]
    U99= [x for x in UprofAvgAbs if x >= 0.99*Ue][0]
    U999= [x for x in UprofAvgAbs if x >= 0.999*Ue][0]
    
    delta = y[UeIndex]
    delta95 = (0.95*Ue-UprofAvg[U95Index-1])/(UprofAvg[U95Index]-UprofAvg[U95Index-1])*(y[U95Index]-y[U95Index-1])+y[U95Index-1]
    delta99 = (0.99*Ue-UprofAvg[U99Index-1])/(UprofAvg[U99Index]-UprofAvg[U99Index-1])*(y[U99Index]-y[U99Index-1])+y[U99Index-1]
    delta999 = (0.999*Ue-UprofAvg[U999Index-1])/(UprofAvg[U999Index]-UprofAvg[U999Index-1])*(y[U999Index]-y[U999Index-1])+y[U999Index-1]
    
    print "Ue = " + str(Ue)
    ## delta* (displacement thickness
    deltas = intsum([(1.0-x/Ue) for x in UprofAvg[0:UeIndex]],[x for x in y[0:UeIndex]])
    deltas2 = np.trapz([(1.0-x/Ue) for x in UprofAvg[0:UeIndex]],x=[x for x in y[0:UeIndex]])
    ## theta  (momentoum thickness
    theta = intsum([(x/Ue)*(1.0-x/Ue) for x in UprofAvg[0:UeIndex]],[x for x in y[0:UeIndex]])
    ##shape factor
    H = deltas/theta
    
    ##skin friction based on inlet velocity
    cf1 = 2*nu*wgu/Uinlet**2
    
    ##skin friction based on edge velocity
    cf2 = 2*nu*wgu/Ue**2
    if cf1==0:
        G=1e10
    else:
        G = (H-1)/H*np.sqrt(2.0/cf2)
        
    if np.isnan(G):
        print "Reverese Flow!!!"
        G = (H-1)/H*np.sqrt(-2.0/cf2)
    utau = np.sqrt(nu*wgu)
    if np.isnan(utau):
        utau=np.sqrt(-nu*wgu)
    
    #writing profiles 
    dSArray = np.ones(len(y))*dS
    UeArray = np.ones(len(y))*Ue
    utauArray = np.ones(len(y))*utau
    deltaArray = np.ones(len(y))*delta
    delta95Array = np.ones(len(y))*delta95
    delta99Array = np.ones(len(y))*delta99
    delta999Array = np.ones(len(y))*delta999
    deltasArray = np.ones(len(y))*deltas
    thetaArray = np.ones(len(y))*theta
    pSArray = np.ones(len(y))*pS
    #Temperature
    if tname!=0:
        T1eIndex = len(y) #T1profAvg.index(T1e)
        ## edge temperature
        T1e = T0 # T1profAvg[T1eIndex-1] #T1e = max(T1profAvg)
        
        ## edge temperature
        T2eIndex = len(y) #T2profAvg.index(T1e)
        T2e = T0 #T2profAvg[T2eIndex-1]# T2e = min(T2profAvg)
        #print T2eIndex
        
        DT1 = T1S-T1e
        DT2 = T2S-T2e
    #    print DT1
    #    print DT2
        if DT1==0:
            thetaT1prof = [ 0 for x in T1profAvg]
            DT1=1e-10
        else:
            thetaT1prof = [ (x-T1e)/DT1 for x in T1profAvg]
        if DT2==0:
            thetaT2prof = [ 0 for x in T2profAvg]
            DT2=1e-10
        else:
            thetaT2prof = [(x-T2e)/DT2 for x in T2profAvg]
            
        St1 = alfa * gradT / (Ue*(T1S-T0))
        St2 = alfa * gradT / (Ue*(T1S-T0))
        
        #    print thetaT1prof
        #    print DT1
        #    print DT2
        
        #thetaT1prof = [(x-T1e)/DT1 for x in T1profAvg[0:T1eIndex]]
        #thetaT2prof = [(x-T2e)/DT2 for x in T2profAvg[0:T2eIndex]]
        #    print thetaT1prof
        #print [i for i in range(len(thetaT1prof))][0]
        
        T195Index = [i for i in range(len(thetaT1prof)) if thetaT1prof[i] <= 0.05][0]
        T199Index = [i for i in range(len(thetaT1prof)) if thetaT1prof[i] <= 0.01][0]
        T1999Index = [i for i in range(len(thetaT1prof)) if thetaT1prof[i] <= 0.005][0]
        
        
        T295Index = [i for i in range(len(thetaT2prof)) if thetaT2prof[i] <= 0.05][0]
        T299Index = [i for i in range(len(thetaT2prof)) if thetaT2prof[i] <= 0.01][0]
        T2999Index = [i for i in range(len(thetaT2prof)) if thetaT2prof[i] <= 0.005][0]
        
        
        deltaT1 = y[T1eIndex-1]
        deltaT195 = y[T195Index]
        deltaT199 = y[T199Index]
        deltaT1999 = y[T1999Index]
        
        
        
        deltaT2 = y[T2eIndex-1]
        deltaT295 = y[T295Index]
        deltaT299 = y[T299Index]
        deltaT2999 = y[T2999Index]
        
        
        deltaT1W = intsum([(x-T1e)/DT1 for x in T1profAvg[0:T1eIndex]],[x for x in y[0:T1eIndex]])
        deltaT2W = intsum([(x-T2e)/DT2 for x in T2profAvg[0:T2eIndex]],[x for x in y[0:T2eIndex]])
        
        
        
        T1eArray = np.ones(len(y))*T1e
        T2eArray = np.ones(len(y))*T2e
        
        T1SArray = np.ones(len(y))*T1S
        T2SArray = np.ones(len(y))*T2S
        
        pSArray = np.ones(len(y))*pS
        
        deltaT1WArray = np.ones(len(y))*deltaT1W
        deltaT2WArray = np.ones(len(y))*deltaT2W
        
        delta95T1Array = np.ones(len(y))*deltaT195
        delta95T2Array = np.ones(len(y))*deltaT295
        
        delta99T1Array = np.ones(len(y))*deltaT199
        delta99T2Array = np.ones(len(y))*deltaT299
        
    else:
        T1eIndex = 0 #T1profAvg.index(T1e)
        ## edge temperature
        T1e = T0 # T1profAvg[T1eIndex-1] #T1e = max(T1profAvg)
        
        ## edge temperature
        T2eIndex = 0 #T2profAvg.index(T1e)
        T2e = T0 #T2profAvg[T2eIndex-1]# T2e = min(T2profAvg)
        #print T2eIndex
        
        thetaT1prof = [ 0 for x in T1profAvg]
        thetaT2prof = [ 0 for x in T1profAvg]
        
        St1 = 0
        St2 = 0
        
        #    print thetaT1prof
        #    print DT1
        #    print DT2
        
        #thetaT1prof = [(x-T1e)/DT1 for x in T1profAvg[0:T1eIndex]]
        #thetaT2prof = [(x-T2e)/DT2 for x in T2profAvg[0:T2eIndex]]
        #    print thetaT1prof
        #print [i for i in range(len(thetaT1prof))][0]
        
        
        deltaT1 = 0
        deltaT195 = 0
        deltaT199 = 0
        deltaT1999 = 0
        
        
        
        deltaT2 = 0
        deltaT295 = 0
        deltaT299 = 0
        deltaT2999 = 0
        
        
        deltaT1W = 0
        deltaT2W = 0
        #writing profiles 
                
        T1eArray = np.ones(len(y))*0
        T2eArray = np.ones(len(y))*0
        
        T1SArray = np.ones(len(y))*0
        T2SArray = np.ones(len(y))*0
        

        
        deltaT1WArray = np.ones(len(y))*0
        deltaT2WArray = np.ones(len(y))*0
        
        delta95T1Array = np.ones(len(y))*0
        delta95T2Array = np.ones(len(y))*0
        
        delta99T1Array = np.ones(len(y))*0
        delta99T2Array = np.ones(len(y))*0
    
    #print len(thetaT2prof)
    #print len(T2profAvg)
    
    
    #    print pS
    #    print pSArray
    np.savetxt(Profiles,np.c_[  dSArray,UeArray,utauArray,deltaArray,delta95Array,delta99Array,delta999Array,deltasArray,thetaArray, \
                                T1eArray,T2eArray,T1SArray,T2SArray,\
                                pSArray, \
                                deltaT1WArray,deltaT2WArray,\
                                delta95T1Array,delta95T2Array,\
                                delta99T1Array,delta99T2Array,\
                                y,\
                                pprofAvg,\
                                UprofAvg,VprofAvg, WprofAvg,\
                                UUprofAvg,VVprofAvg,WWprofAvg,UVprofAvg,UWprofAvg,VWprofAvg,\
                                gradUprofAvg,\
                                nuSgsprofAvg,\
                                T1profAvg,T2profAvg,\
                                thetaT1prof, thetaT2prof,\
                                TT1profAvg,TT2profAvg,\
                                VT1profAvg,VT2profAvg,\
                                UT1profAvg,UT2profAvg,\
                                TuprofAvg,\
                                gradT1profAvg,gradT2profAvg],
                            fmt="%15.12f "*46)
    
    
    print "delta  = " + str(delta)
    print "deltas = " + str(deltas)
    print "deltaT1W = " + str(deltaT1W)
    print "deltaT199 = " + str(deltaT199)
    print "deltaT299 = " + str(deltaT299)
    print "deltaT1W = " + str(deltaT1W)
                                          
    print "theta = " + str(theta)
    print "shape Factor = " +str(deltas/theta)
    
    ### statistics
    #flow
    UUEdge = UUprofAvg[U99Index]
    VVEdge = VVprofAvg[U99Index]
    WWEdge = WWprofAvg[U99Index]
    
    UVEdge = UVprofAvg[U99Index]
    UWEdge = UWprofAvg[U99Index]
    VWEdge = VWprofAvg[U99Index]
    
    TuEdge = TuprofAvg[U99Index]
    
    UUMax = max(UUprofAvg)
    VVMax = max(VVprofAvg)
    WWMax = max(WWprofAvg)
    
    UVMax = min(UVprofAvg)
    UWMax = max(UWprofAvg)
    VWMax = max(VWprofAvg)
    
    TuMax = max(TuprofAvg)
    
    UUMaxIndex = UUprofAvg.index(UUMax)
    VVMaxIndex = VVprofAvg.index(VVMax)
    WWMaxIndex = WWprofAvg.index(WWMax)
    
    UVMaxIndex = UVprofAvg.index(UVMax)
    UWMaxIndex = UWprofAvg.index(UWMax)
    VWMaxIndex = VWprofAvg.index(VWMax)
    
    TuMaxIndex = TuprofAvg.index(TuMax)
    
    yUUMax = y[UUMaxIndex]
    yVVMax = y[VVMaxIndex]
    yWWMax = y[WWMaxIndex]
    yUVMax = y[UVMaxIndex]
    yUWMax = y[UWMaxIndex]
    yVWMax = y[VWMaxIndex]
    
    yTuMax = y[TuMaxIndex]
    
    UVMaxNon = UVMax / Ue**2
    #thermal
    TT1Max = max(TT1profAvg)
    TT2Max = max(TT2profAvg)
    
    TT1MaxIndex = TT1profAvg.index(TT1Max)
    TT2MaxIndex = TT2profAvg.index(TT2Max)
    if TT1Max ==0:
        yTT1Max =0.0
    else:
        yTT1Max = y[TT1MaxIndex]
    
    if TT2Max ==0:
        yTT2Max =0.0
    else:
        yTT2Max = y[TT2MaxIndex]
    
    #print "UVIndex = " + str(UVMaxIndex)
    #print "TuxIndex = " + str(UUMaxIndex)
    #print "UVmax = " + str(UVMax)
    #print "UUMax = " + str(UUMax)
    #print "y @UVMAx = " +str(yUVMax)
    #print "y @UUMax = " +str(yUUMax)
    
    #print [x for x in UprofAvg if x > 0.99*Ue]
    #print [i for i in range(len(UprofAvg)) if UprofAvg[i] > 0.99*Ue]
    #print y[0:UeIndex+1]
    
    
    boundarylayer.write( ("%15.12f "*34 +" \n") % \
                    (dS, Ue, pS, St1, St2, T1S, T2S, T1e, T2e, delta, \
                    delta95, delta99, delta999, deltas, theta, \
                    H, G, \
                    deltaT1, deltaT195, deltaT199, deltaT1999, \
                    deltaT2, deltaT295, deltaT299, deltaT2999, \
                    deltaT1W, deltaT2W, cf1, cf2, utau, ri, UUEdge, VVEdge, WWEdge  ))
    
    str1="%15.12f "*29 + " \n"
    #print str1
    Maxvalues.write( ("%15.12f "*36 +" \n") % \
            (dS, Ue, T1S, T2S, T1e,\
            T2e, UUMax, VVMax, WWMax, UVMax,\
            UWMax, VWMax, UVMaxNon, TuMax, UUEdge, VVEdge,\
            WWEdge, UVEdge, UWEdge, VWEdge, TuEdge, TT1Max,\
            TT2Max,yUUMax, yVVMax, yWWMax, yUVMax,\
            yUWMax,yVWMax, yTuMax, yTT1Max, yTT2Max, delta, delta95, delta99, delta999))
    
    index=index+1
    #print "index=" + str(index)
    del y
    del UprofAvg
    del UprofAvgAbs
    del VprofAvg
    del WprofAvg
    del T1profAvg
    del T2profAvg
    del TT1profAvg
    del TT2profAvg
    #    del UT1profAvg
    #    del UT2profAvg
    del UUprofAvg
    del VVprofAvg
    del WWprofAvg
    del UVprofAvg
    del UWprofAvg
    del VWprofAvg
    del TuprofAvg
    del dSArray
    del UeArray
    del utauArray
    
    del T1eArray
    del T2eArray
    
    del T1SArray
    del T2SArray
boundarylayer.close()
Maxvalues.close()
Profiles.close()
del Profiles



Delete(integrateVariables2)
del integrateVariables2
Delete(slice2)
del slice2
Delete(clip1)
del clip1
Delete(clip2)
del clip2

Delete(integrateVariables1)
del integrateVariables1
Delete(plotOverLine1)
del plotOverLine1 
                                                                                                                                                                                                                                                                                                        
