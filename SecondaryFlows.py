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
mm = 400
first=1e-5
maxy=0.01
l=np.logspace(np.log10(first), np.log10(maxy), num=mm,base=10)
l=np.insert(l,0,0)
#print l
nn = 2
Smin=0
Smax=0.6
deltaS=0.01
nS=(Smax-Smin)/deltaS+1
print nS
SS=np.linspace(Smin,Smax,nS)
print SS
print len(SS)
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

cdata= np.loadtxt('D:/Yousef/paraview-scripts/curvecascade')
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
caseName = sourceName # raw_input("Please enter case name: ")
source = FindSource(sourceName + '.foam')
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
slice2 = Slice(Input=source)
slice2.SliceType = 'Plane'
slice2.SliceOffsetValues = [0.0]
slice2.SliceType.Origin = [0, 0, 0.18]
slice2.SliceType.Normal = [0.0,0.0, 1.0]


# create a new 'Slice'
slice1 = Slice(Input=source)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [xi, yi, 0]
slice1.SliceType.Normal = [X2-X1,Y2-Y1, 0.0]



surfaceNormals = GenerateSurfaceNormals(Input=slice1)
# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

clip1 = Clip(Input=surfaceNormals)
clip1.ClipType = 'Plane'
clip1.ClipType.Origin = [clipx, clipy, 0.127]
clip1.ClipType.Normal = [-(Y2-Y1),(X2-X1),0]
clip1.Invert=0


d = 0.1

r=((Y2-Y1)**2+(X2-X1)**2)**0.5







############################################
##############################
##############################
##############################

####### initializing files #######
path1='I:/PostProcess/pvout/'
path2=".csv"

path = path1 + caseName + "_base" + str(dS) + "_"  + path2
# save data
#SaveData(path, proxy=slice2)
clip1 = GetActiveSource()
index=1
for dS in SS:
    
       
    Sa =  St + dS*sh  #absolute S
    path = path1 + caseName + "_S" + str(dS) + path2
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

    #    if Sa>0:
    clipy=yi-0.001*(X2-X1)/r
    clipx=xi+0.001*(Y2-Y1)/r
    #    else:
    #        clipy=yi+0.01*(X2-X1)/r
    #        clipx=xi-0.01*(Y2-Y1)/r
    
    clip1.ClipType.Origin = [clipx, clipy, 0]
    clip1.ClipType.Normal = [-(Y2-Y1),(X2-X1),0]
    # save data
    SaveData(path, proxy=clip1)
    # get active source.
    
    