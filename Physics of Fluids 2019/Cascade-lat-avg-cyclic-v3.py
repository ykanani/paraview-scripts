##v2 added back the slices to reduce computational efforts
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
mm = 200 
first=1e-5
maxy=0.02
l=np.logspace(np.log10(first), np.log10(maxy), num=mm,base=10)
l=np.insert(l,0,0)
#print l
nn = 2
Smin=0#-0.45
Smax=0.6
deltaS=0.01
nS=(Smax-Smin)/deltaS+1
SS=np.linspace(Smin,Smax,nS)
print SS
L=0.4976

###############################################
n= 50
delta=0.02
zmin=-0.06#0.117 #0.1
zmax=0.06#0.137 #0.154
sp=zmax-zmin  #spanwise distance

###FOR 2D case
#zmin=0.01
#zmax=0.02
Pr = 0.711
nu = 1.55/100000.0
alfa = nu/Pr
gradT = 38000.0
T0=20
Uexit = 15.95
Uinlet = 4
#############
###############################################
os.chdir(os.path.dirname(os.path.realpath(__file__)))
print(os.getcwd())
cdata= np.loadtxt('curvecascade')
cx = cdata[:,0] #curve x coordinates
cy = cdata[:,1] #curve y coordinates
cS = cdata[:,2] #curve S
cm = cdata[:,3] #curve slope
cR = cdata[:,4] #curve slope

sh = 1 # length scale
St = 0 #S leading edge

dS =0.2
Sa =  St + dS*sh  #absolute S

#print cS.size
w = np.where(cS<Sa)
#print w
wtemp = w[0]
S1 = cS[w][wtemp.size-1] #lower bound
S2 = cS[np.where(cS>Sa)][0] #upper bound

X1 = cx[w][wtemp.size-1] #lower bound
X2 = cx[np.where(cS>Sa)][0] #upper bound
Y1 = cy[w][wtemp.size-1] #lower bound
Y2 = cy[np.where(cS>Sa)][0] #upper bound
M1 = cm[w][wtemp.size-1] #lower bound
M2 = cm[np.where(cS>Sa)][0] #upper bound
R1 = cR[w][wtemp.size-1] #lower bound
R2 = cR[np.where(cS>Sa)][0] #upper bound
xi = X1 + (Sa - S1)*(X2-X1)/(S2-S1)
yi = Y1 + (Sa - S1)*(Y2-Y1)/(S2-S1)
mi = M1 + (Sa - S1)*(M2-M1)/(S2-S1)
ri = R1 + (Sa - S1)*(R2-R1)/(S2-S1)
r=((Y2-Y1)**2+(X2-X1)**2)**0.5

##############################
##############################
############# 1.55/100000/0.711*38000/(15.8*(T1Mean-20))


###############################
## find source
sourceName = raw_input("Please enter source name: ")
tname =input("Please enter 0 for T, 1 for T1, and 2 for T1 and T2")
#sourceName = "CL2"
#tname =input("Please enter 0 for T1, 1 for T2")
caseName = sourceName # raw_input("Please enter case name: ")
print sourceName
source = FindSource(sourceName + '.foam')
sourceS = FindSource(sourceName + 'S.foam')
print source
# get active view
renderView1 = GetActiveViewOrCreate('RenderView')


# init the 'Plane' selected for 'SliceType'
if Sa>0:
    clipy=yi-0.001*(X2-X1)/r
    clipx=xi+0.001*(Y2-Y1)/r
else:
    clipy=yi+0.001*(X2-X1)/r
    clipx=xi-0.001*(Y2-Y1)/r

# create a new 'Slice'
slice1 = Slice(Input=source)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [xi, yi, 0]
slice1.SliceType.Normal = [X2-X1,Y2-Y1, 0.0]




# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

clip1 = Clip(Input=slice1)
clip1.ClipType = 'Plane'
clip1.ClipType.Origin = [clipx, clipy, 0.127]
clip1.ClipType.Normal = [-(Y2-Y1),(X2-X1),0]
clip1.Invert=0

plotOverLine1 = PlotOverLine(Input=clip1,
    Source='High Resolution Line Source')
#start = input("Please enter start: ")
#end = input("Please enter end: ")
#sp=0.031750
# init the 'High Resolution Line Source' selected for 'Source'
d = 0.1

r=((Y2-Y1)**2+(X2-X1)**2)**0.5


plotOverLine1.Source.Point1 = [xi-d*(Y2-Y1)/r, yi+d*(X2-X1)/r, zmin]
plotOverLine1.Source.Point2 = [xi-d*(Y2-Y1)/r, yi+d*(X2-X1)/r, zmax]

# Properties modified on plotOverLine1
plotOverLine1.Tolerance = 2.22044604925031e-16

# create a new 'Integrate Variables'
integrateVariables1 = IntegrateVariables(Input=plotOverLine1)

#####################################surface

# create a new 'Slice' on the surface

slice2 = Slice(Input=sourceS)
slice2.SliceType = 'Plane'
slice2.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice2.SliceType.Origin = [xi, yi, 0]
slice2.SliceType.Normal = [X2-X1,Y2-Y1, 0.0]



clip4 = Clip(Input=slice2)
clip4.ClipType = 'Plane'

clip4.ClipType.Origin = [clipx, clipy, 0.127]
clip4.ClipType.Normal = [-(Y2-Y1),(X2-X1),0]
clip4.Invert=0

integrateVariables2 = IntegrateVariables(Input=clip4)
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


#SS=[0.01 ,.1 ,.2 ,.22 ]
#SS=[0.1]
#SS=[-10]

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



index=1
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
    
    Sa =  St + dS*sh  #absolute S
    print "Current Location + " + str(Sa)
    w =np.where(cS<Sa)
    wtemp = w[0]
    #print cS
    S1 = cS[w][wtemp.size-1] #lower bound
    S2 = cS[np.where(cS>Sa)][0] #upper bound
    
    X1 = cx[w][wtemp.size-1] #lower bound
    X2 = cx[np.where(cS>Sa)][0] #upper bound
    
    Y1 = cy[w][wtemp.size-1] #lower bound
    Y2 = cy[np.where(cS>Sa)][0] #upper bound
    
    M1 = cm[w][wtemp.size-1] #lower bound
    M2 = cm[np.where(cS>Sa)][0] #upper bound
    
    R1 = cR[w][wtemp.size-1] #lower bound
    R2 = cR[np.where(cS>Sa)][0] #upper bound
    
    
    xi = X1 + (Sa - S1)*(X2-X1)/(S2-S1)
    yi = Y1 + (Sa - S1)*(Y2-Y1)/(S2-S1)
    
    mi = M1 + (Sa - S1)*(M2-M1)/(S2-S1)
    ri = R1 + (Sa - S1)*(R2-R1)/(S2-S1)
    
    r=((Y2-Y1)**2+(X2-X1)**2)**0.5

    # modify slice1
    slice1.SliceType.Origin = [xi, yi, 0]
    slice1.SliceType.Normal = [X2-X1,Y2-Y1, 0.0]
    # modify slice 2
    slice2.SliceType.Origin = [xi, yi, 0]
    slice2.SliceType.Normal = [X2-X1,Y2-Y1, 0.0]
    
    #    if Sa>0:
    clipy=yi-0.001*(X2-X1)/r
    clipx=xi+0.001*(Y2-Y1)/r
    #    else:
    #        clipy=yi+0.01*(X2-X1)/r
    #        clipx=xi-0.01*(Y2-Y1)/r
    
    clip1.ClipType.Origin = [clipx, clipy, 0]
    clip1.ClipType.Normal = [-(Y2-Y1),(X2-X1),0]
    
    clip4.ClipType.Origin = [clipx, clipy, 0]
    clip4.ClipType.Normal = [-(Y2-Y1),(X2-X1),0]

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
    Profiles.write( "ZONE T=\" " + sourceName + "_" + str(index) + "_" +  str(round(dS/L,2)) + "\" \n" )
    j=0
    
    for d in l:
        
        #d = maxy/(mm-1)*j
        #d=1e-6*(1-1.01**j)/(1-1.01)
        #print j
        print "d = " +  str(d)
        plotOverLine1.Source.Point1 = [xi-d*(Y2-Y1)/r, yi+d*(X2-X1)/r, zmin]
        plotOverLine1.Source.Point2 = [xi-d*(Y2-Y1)/r, yi+d*(X2-X1)/r, zmax]
        
        pdi = servermanager.Fetch(integrateVariables1)
        pdiS = servermanager.Fetch(integrateVariables2)
        print " HERE 1 "
        #print pdi
        if j ==0 :
            try:
                wgum = pdiS.GetPointData().GetArray("wallGradUMean")
                pMeanS = pdiS.GetPointData().GetArray("pMean")
                
                if tname==0:
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
                print "No wall!!! skipping this dS/L =" + str(dS/L)
                wgu = -1
                ps = 0
                if tname==0:
                    T1MeanS = pdi.GetPointData().GetArray("TMean")
                    T2MeanS = pdi.GetPointData().GetArray("TMean")
                    T1Prime2MeanS = pdi.GetPointData().GetArray("TPrime2Mean")
                    T2Prime2MeanS = pdi.GetPointData().GetArray("TPrime2Mean")
                    gradT1MeanS = pdi.GetPointData().GetArray("gradTMeantn")
                    gradT2MeanS = pdi.GetPointData().GetArray("gradTMeantn")
                elif tname==1:
                    T1MeanS = pdi.GetPointData().GetArray("T1Mean")
                    T2MeanS = pdi.GetPointData().GetArray("T2Mean")
                    T1Prime2MeanS = pdi.GetPointData().GetArray("T1Prime2Mean")
                    T2Prime2MeanS = pdi.GetPointData().GetArray("T2Prime2Mean")
                    gradT1MeanS = pdi.GetPointData().GetArray("gradT1Meantn")
                    gradT2MeanS = pdi.GetPointData().GetArray("gradT1Meantn")
                
                T1S = T1MeanS.GetTuple(0)[0]/sp
                T2S = T2MeanS.GetTuple(0)[0]/sp
                TT1S = T1Prime2MeanS.GetTuple(0)[0]/sp
                TT2S = T1Prime2MeanS.GetTuple(0)[0]/sp
        
        #reading variables
        pMean = pdi.GetPointData().GetArray("pMean")
        UMean = pdi.GetPointData().GetArray("UMeantn")
        Ums = pdi.GetPointData().GetArray("UPrime2Meantn")
        gradUMean = pdi.GetPointData().GetArray("gradUMeantn")   
        nuSgs = pdi.GetPointData().GetArray("nut")   ##############################FOR LAMINAR I CHANGED NUT TO PMEAN !!!!
        if tname==0:
            T1Mean = pdi.GetPointData().GetArray("TMean")
            T2Mean = pdi.GetPointData().GetArray("TMean")
            
            T1Prime2Mean = pdi.GetPointData().GetArray("TPrime2Mean")
            T2Prime2Mean = pdi.GetPointData().GetArray("TPrime2Mean")
            
            UT1Mean = pdi.GetPointData().GetArray("utAvgtn")
            UT2Mean = pdi.GetPointData().GetArray("utAvgtn")
            
            gradT1Mean = pdi.GetPointData().GetArray("gradTMeantn")
            gradT2Mean = pdi.GetPointData().GetArray("gradTMeantn")
        elif tname==1:
            T1Mean = pdi.GetPointData().GetArray("T1Mean")
            T2Mean = pdi.GetPointData().GetArray("T1Mean")
            
            T1Prime2Mean = pdi.GetPointData().GetArray("T1Prime2Mean")
            T2Prime2Mean = pdi.GetPointData().GetArray("T1Prime2Mean")
            
            UT1Mean = pdi.GetPointData().GetArray("ut1Avgtn")
            UT2Mean = pdi.GetPointData().GetArray("ut1Avgtn")
            
            gradT1Mean = pdi.GetPointData().GetArray("gradT1Meantn")
            gradT2Mean = pdi.GetPointData().GetArray("gradT1Meantn")
        elif tname==2:
            T1Mean = pdi.GetPointData().GetArray("T1Mean")
            T2Mean = pdi.GetPointData().GetArray("T2Mean")
            
            T1Prime2Mean = pdi.GetPointData().GetArray("T1Prime2Mean")
            T2Prime2Mean = pdi.GetPointData().GetArray("T2Prime2Mean")
            
            UT1Mean = pdi.GetPointData().GetArray("ut1Avgtn")
            UT2Mean = pdi.GetPointData().GetArray("ut2Avgtn")
            
            gradT1Mean = pdi.GetPointData().GetArray("gradT1Meantn")
            gradT2Mean = pdi.GetPointData().GetArray("gradT2Meantn")    
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
            
            gradUprofAvg.append(gradUMean.GetTuple(0)[3]/sp)
            
            nuSgsprofAvg.append(nuSgs.GetTuple(0)[0]/sp)
            
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
    delta95 = y[U95Index]
    delta99 = y[U99Index]
    delta999 = y[U999Index]
    
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
    
    #Temperature
    
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
        
    St1 = alfa * gradT / (Uexit*(T1S-T0))
    St2 = alfa * gradT / (Uexit*(T1S-T0))
    
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
Delete(clip4)
del clip4


Delete(integrateVariables1)
del integrateVariables1
Delete(plotOverLine1)
del plotOverLine1 
                                                                                                                                                                                                                                                                                                        
