##v2 added back the slices to reduce computational efforts
##v3 added other components of gradT
##v4 converted lists to numpy array
#import sys
#sys.modules[__name__].__dict__.clear()
#### import the simple module from the paraview
from paraview.simple import *
import numpy as np
import os, sys

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
Smin=-0.45
Smax=0.6
deltaS=0.01
nS=(Smax-Smin)/deltaS+1
SS=np.linspace(Smin,Smax,nS)
print SS
L=0.4976

###############################################
n= 50
delta=0.02
zmin=-0.05999#0.117 #0.1
zmax=0.05999#0.137 #0.154
sp=zmax-zmin  #spanwise distance
print("sp=" + str(sp))
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

location =4 #input("Please enter the location that you are running this script, 0:Marvericks,1:Office PC,2:BW,3:Stampede2 4:Neso :")
if location==0:
	path1="/home/03624/ykanani/tempparaview/"
elif location==1:
	path1='E:/PostProcess/pvout/'
elif location==2:
	path1='/u/sciteam/kanani/tempparaview/'
elif location==3:
	path1='/home1/03624/ykanani/tempParaview/'
elif location==4:
    path1='I:/PostProcess/pvout/'


os.chdir(path1)
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
sourceName = 'ccc9' #raw_input("Please enter source name: ")
tname =2 #input("Please enter 0 for T, 1 for T1, and 2 for T1 and T2")
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

path2=".plt"
print(path1)

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

def intsum(f,x):
        #print "in intsum ..........."
        #print x
        #print dx
        lenx = len(x)
        #print lenx
        f0 = f[0:lenx-1]
        f1 = f[1:lenx]
        fAve = (f0+f1)/2 #[sum(x)/2.0 for x in zip(x0,x1)]
        
        dx = np.diff(x) ##[xx1 - xx0 for xx1 , xx0 in zip(dx1,dx0)]
        
        return np.inner(fAve,dx)



index=1
for dS in SS:
    
    # restarting main profile variables
    y=np.zeros(mm)
    pprofAvg=np.zeros(mm)
    
    UprofAvg=np.zeros(mm)
    VprofAvg=np.zeros(mm)
    WprofAvg=np.zeros(mm)
    
    UVprofAvg=np.zeros(mm)
    VWprofAvg=np.zeros(mm)
    UWprofAvg=np.zeros(mm)
    
    UUprofAvg=np.zeros(mm)
    VVprofAvg=np.zeros(mm)
    WWprofAvg=np.zeros(mm)
    
    gradUprofAvg=np.zeros(mm)
    
    nuSgsprofAvg=np.zeros(mm)
    
    T1profAvg=np.zeros(mm)
    T2profAvg=np.zeros(mm)
    
    TT1profAvg=np.zeros(mm)
    TT2profAvg=np.zeros(mm)
    
    VT1profAvg=np.zeros(mm)
    VT2profAvg=np.zeros(mm)

    UT1profAvg=np.zeros(mm)
    UT2profAvg=np.zeros(mm)

    TuprofAvg=np.zeros(mm)

    gradT1NprofAvg=np.zeros(mm)
    gradT2NprofAvg=np.zeros(mm)
    gradT1TprofAvg=np.zeros(mm)
    gradT2TprofAvg=np.zeros(mm)

    
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
                                gradT1n gradT2n \
                                gradT1t gradT2t \
                                gradT1nS gradT2nS \
                                gradT1tS gradT2tS \
                                " + " \n" )
    Profiles.write( "ZONE T=\" " + sourceName + "_" + str(index) + "_" +  str(round(dS/L,2)) + "\" \n" )
    #j=0
    #for d in l:
    for j in range(mm):
        d=l[j]
        #d = maxy/(mm-1)*j
        #d=1e-6*(1-1.01**j)/(1-1.01)
        #print j
        #print "d = " +  str(d)
        plotOverLine1.Source.Point1 = [xi-d*(Y2-Y1)/r, yi+d*(X2-X1)/r, zmin]
        plotOverLine1.Source.Point2 = [xi-d*(Y2-Y1)/r, yi+d*(X2-X1)/r, zmax]
        
        pdi = servermanager.Fetch(integrateVariables1)
        pdiS = servermanager.Fetch(integrateVariables2)
        #print " HERE 1 "
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
                gradT1NS = gradT1MeanS.GetTuple(0)[1]/sp
                gradT2NS = gradT2MeanS.GetTuple(0)[1]/sp
                gradT1TS = gradT1MeanS.GetTuple(0)[0]/sp
                gradT2TS = gradT2MeanS.GetTuple(0)[0]/sp
                print("GradT1S0=",gradT1MeanS.GetTuple(0)[0]/sp)
                print("GradT1S1=",gradT1MeanS.GetTuple(0)[1]/sp)
                print("GradT2S0=",gradT2MeanS.GetTuple(0)[0]/sp)
                print("GradT2S1=",gradT2MeanS.GetTuple(0)[1]/sp)
                
                print "WGU =" + str(wgu)
            except:
                print "Unexpected error:", sys.exc_info()[0]
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
        nuSgs = pdi.GetPointData().GetArray("pMean")   ##############################FOR LAMINAR I CHANGED NUT TO PMEAN !!!!
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
            y[j] = 0.0
            pprofAvg[j]=pS
            
            UprofAvg[j]=0.0
            VprofAvg[j]=0.0
            WprofAvg[j]=0.0
            
            UVprofAvg[j]=0.0
            UWprofAvg[j]=0.0
            VWprofAvg[j]=0.0
            
            UUprofAvg[j]=0.0
            VVprofAvg[j]=0.0
            WWprofAvg[j]=0.0
            
            gradUprofAvg[j]=wgu
            
            nuSgsprofAvg[j]=0.0
            
            T1profAvg[j]= T1S
            T2profAvg[j]= T2S
            
            TT1profAvg[j]=TT1S
            TT2profAvg[j]=TT2S
            
            VT1profAvg[j]=0.0
            VT2profAvg[j]=0.0
            
            UT1profAvg[j]=0.0
            UT2profAvg[j]=0.0
            
            TuprofAvg[j]=0.0
            
            gradT1NprofAvg[j]=-38000
            gradT2NprofAvg[j]=gradT2NS
            gradT1TprofAvg[j]=gradT1TS
            gradT2TprofAvg[j]=gradT2TS
        elif  np.isfinite(UMean.GetTuple(0)[0]/sp):
            
            #print "d="+str(d)
            y[j]=d
            #print j
            pprofAvg[j]=pMean.GetTuple(0)[0]/sp
            
            UprofAvg[j]=UMean.GetTuple(0)[0]/sp
            VprofAvg[j]=UMean.GetTuple(0)[1]/sp
            WprofAvg[j]=UMean.GetTuple(0)[2]/sp
            
            UVprofAvg[j]=Ums.GetTuple(0)[3]/sp
            VWprofAvg[j]=Ums.GetTuple(0)[4]/sp
            UWprofAvg[j]=Ums.GetTuple(0)[5]/sp
            
            UUprofAvg[j]=Ums.GetTuple(0)[0]/sp
            VVprofAvg[j]=Ums.GetTuple(0)[1]/sp
            WWprofAvg[j]=Ums.GetTuple(0)[2]/sp
            
            gradUprofAvg[j]=gradUMean.GetTuple(0)[3]/sp
            
            nuSgsprofAvg[j]=nuSgs.GetTuple(0)[0]/sp
            
            T1profAvg[j]=T1Mean.GetTuple(0)[0]/sp
            T2profAvg[j]=T2Mean.GetTuple(0)[0]/sp
            
            TT1profAvg[j]=T1Prime2Mean.GetTuple(0)[0]/sp
            TT2profAvg[j]=T2Prime2Mean.GetTuple(0)[0]/sp
            
            VT1profAvg[j]=UT1Mean.GetTuple(0)[1]/sp
            VT2profAvg[j]=UT2Mean.GetTuple(0)[1]/sp
            
            UT1profAvg[j]=UT1Mean.GetTuple(0)[0]/sp
            UT2profAvg[j]=UT2Mean.GetTuple(0)[0]/sp
            
            TuprofAvg[j]=np.sqrt((Ums.GetTuple(0)[0]/sp+Ums.GetTuple(0)[1]/sp+Ums.GetTuple(0)[2]/sp)/3)
            
            gradT1NprofAvg[j]=gradT1Mean.GetTuple(0)[1]/sp
            gradT2NprofAvg[j]=gradT2Mean.GetTuple(0)[1]/sp

            gradT1TprofAvg[j]=gradT1Mean.GetTuple(0)[0]/sp
            gradT2TprofAvg[j]=gradT2Mean.GetTuple(0)[0]/sp
        else:
            print('Nan Detected')

    #print T1profAvg
    
    
    UprofAvgAbs = np.absolute(UprofAvg)
    Ue = np.amax(UprofAvg)

    
    U95= UprofAvg[np.where(UprofAvg>=0.95*Ue)]
    U99= UprofAvg[np.where(UprofAvg>=0.99*Ue)]
    U999= UprofAvg[np.where(UprofAvg>=0.999*Ue)]

    
    delta = y[np.where(UprofAvg==Ue)]
    delta95 = y[np.where(UprofAvg==0.95*Ue)]
    delta99 = y[np.where(UprofAvg==0.99*Ue)]
    delta999 = y[np.where(UprofAvg==0.999*Ue)]
    
    print "Ue = " + str(Ue)
    ## delta* (displacement thickness
    BLindices = np.where(UprofAvg <= Ue)
    
    deltas = np.trapz(1.0-UprofAvg[BLindices]/Ue,y[BLindices])
    ## theta  (momentoum thickness
    print('delta=' + str(delta))

    theta = np.trapz(UprofAvg[BLindices]/Ue*(1.0-UprofAvg[BLindices]/Ue),y[BLindices])
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
        thetaT1prof = np.zeros(mm)
        DT1=1e-10
    else:
        thetaT1prof = (T1profAvg-T1e)/DT1
    if DT2==0:
        thetaT2prof = np.zeros(mm)
        DT2=1e-10
    else:
        thetaT2prof = (T2profAvg-T2e)/DT2
        
    St1 = alfa * gradT / (Uexit*(T1S-T0))
    St2 = - alfa * gradT2NS / (Uexit*(T2S-T0))
    
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
    deltaT195 =  y[np.where(thetaT1prof<=0.05)][0]
    deltaT199 = y[np.where(thetaT1prof<=0.01)][0]
    deltaT1999 = y[np.where(thetaT1prof<=0.01)][0]
    
    
    
    deltaT2 = y[T2eIndex-1]
    deltaT295 =  y[np.where(thetaT2prof<=0.05)][0]
    deltaT299 = y[np.where(thetaT2prof<=0.01)][0]
    deltaT2999 = y[np.where(thetaT2prof<=0.01)][0]
    
    
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
    
    gradT1NSArray = np.ones(len(y))*gradT1NS
    gradT2NSArray = np.ones(len(y))*gradT2NS
    gradT1TSArray = np.ones(len(y))*gradT1TS
    gradT2TSArray = np.ones(len(y))*gradT2TS
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
                                gradT1NprofAvg,gradT2NprofAvg,\
                                gradT1TprofAvg,gradT2TprofAvg,\
                                gradT1NSArray,gradT2NSArray,\
                                gradT1TSArray,gradT2TSArray    ],
                            fmt="%15.12f "*52)

    
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
                                                                                                                                                                                                                                                                                                        
