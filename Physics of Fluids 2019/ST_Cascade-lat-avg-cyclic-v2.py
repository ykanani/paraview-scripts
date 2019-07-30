## v2 - properly normalized pressure
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
maxy=0.008#0.02
l=np.logspace(np.log10(first), np.log10(maxy), num=mm,base=10)
print(l)
l=np.insert(l,0,0)
#print l
nn = 2
Smin=-0.45
Smax=0.6
deltaS=0.005
nS=(Smax-Smin)/deltaS+1
SS=np.linspace(Smin,Smax,nS)
print SS
L=0.4976 #chord
###############################################
n= 50
delta=0.02
#### to use for endwall cases  first clip to zmin and zmax
zmin=-0.06 # 0.1##0.117 #0.1
zmax=0.06 # 0.154#0.06#0.137 #0.154
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

cdata= np.loadtxt('E:/PostProcess/PythonScripts/curvecascade')
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

# init the 'Plane' selected for 'SliceType'
if Sa>0:
    clipy=yi-0.001*(X2-X1)/r
    clipx=xi+0.001*(Y2-Y1)/r
else:
    clipy=yi+0.001*(X2-X1)/r
    clipx=xi-0.001*(Y2-Y1)/r
###############################
## find source
sourceName = raw_input("Please enter source name: ")
tname =input("Please enter 0 for T, 1 for T1 and T2")
#sourceName = "CL2"
#tname =input("Please enter 0 for T1, 1 for T2")
caseName = sourceName # raw_input("Please enter case name: ")
print sourceName
sourceS = FindSource(sourceName + 'S.foam')

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')


slice2 = Slice(Input=sourceS)
slice2.SliceType = 'Plane'
slice2.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice2.SliceType.Origin = [xi, yi, 0]
slice2.SliceType.Normal = [X2-X1,Y2-Y1, 0.0]



clip4 = Clip(Input=slice2)
clip4.ClipType = 'Plane'

clip4.ClipType.Origin = [clipx, clipy, 0]
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


path = path1 + caseName + "_St"  + path2
boundarylayer = open(path,'a')
boundarylayer.write( "VARIABLES= dS dS/L St1 St2 T1S T2S cf1 utau pS pnon" + " \n")
boundarylayer.write( "ZONE T=\" " + sourceName + "_St" + "\" \n" )



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
dSArray = [] 
dSLArray = []
St1Array = []
St2Array = []
T1SArray = []
T2SArray = []
cf1Array = []
utauArray = []
pSArray = []
pnonArray = []
for dS in SS:
    
    # restarting main profile variables
    y=[]
    
    
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
    
    # modify slice 2
    slice2.SliceType.Origin = [xi, yi, 0]
    slice2.SliceType.Normal = [X2-X1,Y2-Y1, 0.0]
    
    #    if Sa>0:
    clipy=yi-0.001*(X2-X1)/r
    clipx=xi+0.001*(Y2-Y1)/r
    #    else:
    #        clipy=yi+0.01*(X2-X1)/r
    #        clipx=xi-0.01*(Y2-Y1)/r

    
    clip4.ClipType.Origin = [clipx, clipy, 0]
    clip4.ClipType.Normal = [-(Y2-Y1),(X2-X1),0]
    
    j=0
    
    
        
    pdiS = servermanager.Fetch(integrateVariables2)
        #print pdi

    wgum = pdiS.GetPointData().GetArray("wallGradUMean")
    pMeanS = pdiS.GetPointData().GetArray("pMean") #pMean
                
    if tname==0:
        T1MeanS = pdiS.GetPointData().GetArray("TMean")
        T2MeanS = pdiS.GetPointData().GetArray("TMean")

    elif tname==1:
        T1MeanS = pdiS.GetPointData().GetArray("T1Mean")
        T2MeanS = pdiS.GetPointData().GetArray("T1Mean")  #######################################################################

        
    wgu = np.sqrt(wgum.GetTuple(0)[0]**2+wgum.GetTuple(0)[1]**2)/sp
    pS = pMeanS.GetTuple(0)[0]/sp
    
    T1S = T1MeanS.GetTuple(0)[0]/sp
    T2S = T2MeanS.GetTuple(0)[0]/sp
    #print str(T1S)
    
    print "WGU =" + str(wgu)
            
        
        
    ##skin friction based on inlet velocity
    cf1 = 2*nu*wgu/Uinlet**2
    
    ##skin friction based on edge velocity
    #cf2 = 2*nu*wgu/Ue**2
    
    utau = np.sqrt(nu*wgu)
    if np.isnan(utau):
        utau=np.sqrt(-nu*wgu)
    
    #Temperature
    
    
    St1 = alfa * gradT / (Uexit*(T1S-T0))
    St2 = alfa * gradT / (Uexit*(T1S-T0))
    
    dSArray.append(dS)
    dSLArray.append(dS/L)
    St1Array.append(St1)
    St2Array.append(St2)
    T1SArray.append(T1S)
    T2SArray.append(T2S)
    cf1Array.append(cf1)
    utauArray.append(utau)
    pSArray.append(pS)
    
    
    index=index+1
    #print St1Array
    #print "index=" + str(index)
    
Pmax=max([float(i) for i in pSArray])
pnonArray = [(Pmax-x)/(0.5*Uexit**2) for x in pSArray]   
print len(pnonArray)
print len(pSArray) 
np.savetxt(boundarylayer,np.c_[ dSArray, dSLArray,St1Array,St2Array,T1SArray,T2SArray,cf1Array,utauArray,pSArray,pnonArray],
                            fmt="%15.12f "*10)
boundarylayer.close()





Delete(integrateVariables2)
del integrateVariables2
Delete(slice2)
del slice2
Delete(clip4)
del clip4

                                                                                                                                                                                                                                                                                                        
