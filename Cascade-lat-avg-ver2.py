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
mm = 500
first=1e-7
maxy=0.02
l=np.logspace(np.log10(first), np.log10(maxy), num=mm,base=10)
l=np.insert(l,0,0)
#print l
nn = 2
Smin=-.45
Smax=0.60
deltaS=0.01
nS=(Smax-Smin)/deltaS+1
SS=np.linspace(Smin,Smax,nS)
print SS

###############################################
n= 5
delta=0.02
zmin=0.1
zmax=0.154
sp=zmax-zmin  #spanwise distance

###FOR 2D case
#zmin=0.01
#zmax=0.02
Pr = 0.711
nu = 1.55/100000.0
alfa = nu/Pr
gradT = 38000.0
T0=20
Uexit = 15.6
Uinlet = 4
#############
###############################################

cdata= np.loadtxt('/home/03624/ykanani/paraviewscripts/curvecascade')
cx = cdata[:,0] #curve x coordinates
cy = cdata[:,1] #curve y coordinates
cS = cdata[:,2] #curve S
cm = cdata[:,3] #curve slope


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
xi = X1 + (Sa - S1)*(X2-X1)/(S2-S1)
yi = Y1 + (Sa - S1)*(Y2-Y1)/(S2-S1)
mi = M1 + (Sa - S1)*(M2-M1)/(S2-S1)

r=((Y2-Y1)**2+(X2-X1)**2)**0.5

##############################
##############################
############# 1.55/100000/0.711*38000/(15.8*(T1Mean-20))


###############################
## find source
sourceName = raw_input("Please enter source name: ")
tname =input("Please enter 0 for T, 1 for T1 and T2")
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
slice1.SliceType.Normal = [X2-X1,Y2-Y1, 0.0]




# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

clip1 = Clip(Input=slice1)
clip1.ClipType = 'Plane'

# init the 'Plane' selected for 'SliceType'
if Sa>0:
    clipy=yi-0.001*(X2-X1)/r
    clipx=xi+0.001*(Y2-Y1)/r
else:
    clipy=yi+0.001*(X2-X1)/r
    clipx=xi-0.001*(Y2-Y1)/r

clip1.ClipType.Origin = [clipx, clipy, 0.127]
clip1.ClipType.Normal = [-(Y2-Y1),(X2-X1),0]


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
clip2 = Clip(Input=sourceS)
clip2.ClipType = 'Plane'

# init the 'Plane' selected for 'SliceType'
clip2.ClipType.Origin = [0, 0, zmin]
clip2.ClipType.Normal = [0,0,1]

clip3 = Clip(Input=clip2)
clip3.ClipType = 'Plane'

# init the 'Plane' selected for 'SliceType'
clip3.ClipType.Origin = [0, 0, zmax]
clip3.ClipType.Normal = [0,0,-1]

slice2 = Slice(Input=clip3)
slice2.SliceType = 'Plane'
slice2.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice2.SliceType.Origin = [xi, yi, 0]
slice2.SliceType.Normal = [X2-X1,Y2-Y1, 0.0]



clip4 = Clip(Input=slice2)
clip4.ClipType = 'Plane'

clip4.ClipType.Origin = [clipx, clipy, 0.127]
clip4.ClipType.Normal = [-(Y2-Y1),(X2-X1),0]


integrateVariables2 = IntegrateVariables(Input=clip4)
############################################
##############################
##############################
##############################

####### initializing files #######
path1="/home/03624/ykanani/tempparaview/"
path2=".plt"


path = path1 + caseName + "_BL_All"  + path2
foutv11 = open(path,'a')
foutv11.write( "VARIABLES= dS Ue pS  St1 St2 T1S T2S T1E T2E delta delta95 delta99 delta999 \
                     deltas theta H G \
                     deltaT1 deltaT195 deltaT199 deltaT1999 \
                     deltaT2 deltaT295 deltaT299 deltaT2999 \
                     deltaT1W deltaT2W \
                     cf utau" + " \n")
foutv11.write( "ZONE T=\" " + sourceName + "_BL_ALL" + "\" \n" )

path = path1 + caseName + "_MAX_All"  + path2
foutv12 = open(path,'a')
foutv12.write( "VARIABLES= dS Ue  T1S T2S T1E T2E TuxMax TuyMax TuzMax UVMax UWMax VWMax UVMaxNon,\
TuxEdge, TuyEdge,TuzEdge, UVEdge, UWEdge, VWEdge, \
T1rmsMax T2rmsMax yTuxMax yTuyMax yTuzMax yUVMax yUWMax yVWMax yT1rmsMax yT2rmsMax" + " \n" )
foutv12.write( "ZONE T=\" " + sourceName + "_MAX_ALL" + "\" \n" )

###saving PROFILES
path = path1 + caseName + "_Profiles" + path2
foutv2 = open(path,'a')


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
index=1
for dS in SS:


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

    xi = X1 + (Sa - S1)*(X2-X1)/(S2-S1)
    yi = Y1 + (Sa - S1)*(Y2-Y1)/(S2-S1)

    mi = M1 + (Sa - S1)*(M2-M1)/(S2-S1)

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

    clip1.ClipType.Origin = [clipx, clipy, 0.127]
    clip1.ClipType.Normal = [-(Y2-Y1),(X2-X1),0]
    
    clip4.ClipType.Origin = [clipx, clipy, 0.127]
    clip4.ClipType.Normal = [-(Y2-Y1),(X2-X1),0]


    foutv2.write( "VARIABLES= dS d Ue deltas theta deltaT1W deltaT2W deltaT1 deltaT2 pS  T1E T1S T2E T2S utau U V W T1Mean T2Mean ThetaT1 ThetaT2 \
                            T1Prime2Mean T2Prime2Mean \
                            TuXX TuYY TuZZ UV UW VW" + " \n" )


    foutv2.write( "ZONE T=\" " + sourceName + "_" + str(index) + "\" \n" )
    j=0
    
    for d in l:

                #d = maxy/(mm-1)*j
		#d=1e-6*(1-1.01**j)/(1-1.01)
                #print j
                #print "d = " +  str(d)
                plotOverLine1.Source.Point1 = [xi-d*(Y2-Y1)/r, yi+d*(X2-X1)/r, zmin]
                plotOverLine1.Source.Point2 = [xi-d*(Y2-Y1)/r, yi+d*(X2-X1)/r, zmax]



                pdi = servermanager.Fetch(integrateVariables1)
                pdiS = servermanager.Fetch(integrateVariables2)
                #print pdi
                if j ==0 :
                    try:
                        wgum = pdiS.GetPointData().GetArray("wallGradUMean")
                        pMeanS = pdiS.GetPointData().GetArray("p")

			if tname==0:
                            T1MeanS = pdiS.GetPointData().GetArray("TMean")
                            T2MeanS = pdiS.GetPointData().GetArray("TMean")
                        elif tname==1:
                            T1MeanS = pdiS.GetPointData().GetArray("T1Mean")
                            T2MeanS = pdiS.GetPointData().GetArray("T2Mean")
                        wgu = np.sqrt(wgum.GetTuple(0)[0]**2+wgum.GetTuple(0)[1]**2)/sp
                        pS = pMeanS.GetTuple(0)[0]/sp

			T1S = T1MeanS.GetTuple(0)[0]/sp
                        T2S = T2MeanS.GetTuple(0)[0]/sp
                        print "WGU =" + str(wgu)
                    except:
                        #not on the wall
                        print "No wall!!! skipping this dS =" + str(dS)
                        wgu = -1
                        if tname==0:
                            T1MeanS = pdi.GetPointData().GetArray("TMean")
                            T2MeanS = pdi.GetPointData().GetArray("TMean")
                        elif tname==1:
                            T1MeanS = pdi.GetPointData().GetArray("T1Mean")
                            T2MeanS = pdi.GetPointData().GetArray("T2Mean")
                        pS = p.GetTuple(0)[0]/sp

                        T1S = T1Mean.GetTuple(0)[0]/sp
                        T2S = T2Mean.GetTuple(0)[0]/sp
                #reading variables
                pMean = pdi.GetPointData().GetArray("pMean")

		UMean = pdi.GetPointData().GetArray("UMeantn")
                Ums = pdi.GetPointData().GetArray("UPrime2Meantn")
#               Ept = pdi.GetPointData().GetArray("Ept")
#               upMean = pdi.GetPointData().GetArray("upMeanRot")

                if tname==0:
                    T1Mean = pdi.GetPointData().GetArray("TMean")
                    T1Prime2Mean = pdi.GetPointData().GetArray("TPrime2Mean")
#                   ut1Mean = pdi.GetPointData().GetArray("utMean")
                    T2Mean = pdi.GetPointData().GetArray("TMean")
                    T2Prime2Mean = pdi.GetPointData().GetArray("TPrime2Mean")
#                   ut2Mean = pdi.GetPointData().GetArray("utMean")
                elif tname==1:
                    T1Mean = pdi.GetPointData().GetArray("T1Mean")
                    T1Prime2Mean = pdi.GetPointData().GetArray("T1Prime2Mean")
#                   ut1Mean = pdi.GetPointData().GetArray("ut1MeanRot")

                    T2Mean = pdi.GetPointData().GetArray("T2Mean")
                    T2Prime2Mean = pdi.GetPointData().GetArray("T2Prime2Mean")
#                   ut2Mean = pdi.GetPointData().GetArray("ut2MeanRot")
                #elif tname==2:
                #    TMean = pdi.GetPointData().GetArray("T2Mean")
                #    TPrime2Mean = pdi.GetPointData().GetArray("T2Prime2Mean")
                #    utMean = pdi.GetPointData().GetArray("ut2Mean")


                if d==0:
			#print "if d =0"
                        y = [0.0]
			pprofAvg=[0.0]
                        UprofAvg= [0.0]
                        VprofAvg=[0.0]
                        WprofAvg=[0.0]

                        T1rmsprofAvg=[0.0]
                        T1profAvg= [T1S]
#                        UT1profAvg=[0.0]

                        T2rmsprofAvg=[0.0]
                        T2profAvg= [T2S]
#                        UT2profAvg=[0.0]

                        TUXprofAvg=[0.0]
                        TUYprofAvg=[0.0]
                        TUZprofAvg=[0.0]
                        UVprofAvg=[0.0]
                        UWprofAvg=[0.0]
                        VWprofAvg=[0.0]
                elif  np.isfinite(UMean.GetTuple(0)[0]/sp):
                    try:
			#print "try"
                        y.append(d)
                        pprofAvg.append(pMean.GetTuple(0)[0]/sp)

			UprofAvg.append(UMean.GetTuple(0)[0]/sp)
                        VprofAvg.append(UMean.GetTuple(0)[1]/sp)
                        WprofAvg.append(UMean.GetTuple(0)[2]/sp)

                        T1rmsprofAvg.append(T1Prime2Mean.GetTuple(0)[0]/sp)
                        T1profAvg.append(T1Mean.GetTuple(0)[0]/sp)
#                        UT1profAvg.append(ut1Mean.GetTuple(0)[0]/sp)

                        T2rmsprofAvg.append(T2Prime2Mean.GetTuple(0)[0]/sp)
                        T2profAvg.append(T2Mean.GetTuple(0)[0]/sp)
#                        UT2profAvg.append(ut2Mean.GetTuple(0)[0]/sp)

                        TUXprofAvg.append(np.sqrt(Ums.GetTuple(0)[0]/sp))
                        TUYprofAvg.append(np.sqrt(Ums.GetTuple(0)[1]/sp))
                        TUZprofAvg.append(np.sqrt(Ums.GetTuple(0)[2]/sp))
                        UVprofAvg.append(Ums.GetTuple(0)[3]/sp)
                        VWprofAvg.append(Ums.GetTuple(0)[4]/sp)
                        UWprofAvg.append(Ums.GetTuple(0)[5]/sp)


                    except NameError:
#
                	#print "name error"        
                        y = [d]
                        pprofAvg= [pMean.GetTuple(0)[0]/sp]
                        UprofAvg= [UMean.GetTuple(0)[0]/sp]
                        VprofAvg=[UMean.GetTuple(0)[1]/sp]
                        WprofAvg=[UMean.GetTuple(0)[2]/sp]
                        T1rmsprofAvg=[T1Prime2Mean.GetTuple(0)[0]/sp]
                        T1profAvg= [T1Mean.GetTuple(0)[0]/sp]
#                        UT1profAvg=[ut1Mean.GetTuple(0)[1]/sp]

                        T2rmsprofAvg=[T2Prime2Mean.GetTuple(0)[0]/sp]
                        T2profAvg= [T2Mean.GetTuple(0)[0]/sp]
#                        UT2profAvg=[ut2Mean.GetTuple(0)[1]/sp]

                        TUXprofAvg=[np.sqrt(Ums.GetTuple(0)[0]/sp)]
                        TUYprofAvg=[np.sqrt(Ums.GetTuple(0)[1]/sp)]
                        TUZprofAvg=[np.sqrt(Ums.GetTuple(0)[2]/sp)]
                        UVprofAvg=[Ums.GetTuple(0)[3]/sp]
	                VWprofAvg=[Ums.GetTuple(0)[4]/sp]
			UWprofAvg=[Ums.GetTuple(0)[5]/sp]
		j=j+1
	#	print UprofAvg[-1]
    	#	if UprofAvg[-1]<0:
    	#		whaaaaaat = raw_input("whaaaaaaaaaaaaat ")




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

    ###    velocity
    #print UprofAvg
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

    deltas = intsum([(1.0-x/Ue) for x in UprofAvg[0:UeIndex]],[x for x in y[0:UeIndex]])
    deltas2 = np.trapz([(1.0-x/Ue) for x in UprofAvg[0:UeIndex]],x=[x for x in y[0:UeIndex]])
    theta = intsum([(x/Ue)*(1.0-x/Ue) for x in UprofAvg[0:UeIndex]],[x for x in y[0:UeIndex]])
    H = deltas/theta
    
    cf = 2*nu*wgu/Uinlet**2

    if cf==0:
        G=1e10
    else:
        G = (H-1)/H*np.sqrt(2.0/cf)

    if np.isnan(G):
        print "Reverese Flow!!!"
        G = (H-1)/H*np.sqrt(-2.0/cf)

    utau = np.sqrt(1.62/100000*wgu)
    if np.isnan(utau):
        utau=np.sqrt(-1.62/100000*wgu)
    #writing profiles
    dSArray = np.ones(len(y))*dS
    UeArray = np.ones(len(y))*Ue
    utauArray = np.ones(len(y))*utau
    deltasArray = np.ones(len(y))*deltas
    thetaArray = np.ones(len(y))*theta


    #Temperature

    T1eIndex = len(y) #T1profAvg.index(T1e)
    T1e = T1profAvg[T1eIndex-1] #T1e = max(T1profAvg)



    T2eIndex = len(y) #T2profAvg.index(T1e)
    T2e = T2profAvg[T2eIndex-1]# T2e = min(T2profAvg)
    #print T2eIndex

    DT1 = T1S-T1e
    DT2 = T2S-T2e
#    print DT1
#    print DT2
    if DT1==0:
        thetaT1prof = [ 0 for x in T1profAvg[0:T1eIndex]]
        DT1=1e-10
    else:
        thetaT1prof = [ (x-T1e)/DT1 for x in T1profAvg[0:T1eIndex]]
    if DT2==0:
        thetaT2prof = [ 0 for x in T2profAvg[0:T2eIndex]]
        DT2=1e-10
    else:
        thetaT2prof = [(x-T2e)/DT2 for x in T2profAvg[0:T2eIndex]]

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
    T1999Index = [i for i in range(len(thetaT1prof)) if thetaT1prof[i] <= 0.001][0]


    T295Index = [i for i in range(len(thetaT2prof)) if thetaT2prof[i] <= 0.05][0]
    T299Index = [i for i in range(len(thetaT2prof)) if thetaT2prof[i] <= 0.01][0]
    T2999Index = [i for i in range(len(thetaT2prof)) if thetaT2prof[i] <= 0.001][0]


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
    deltaT1Array = np.ones(len(y))*deltaT199
    deltaT2Array = np.ones(len(y))*deltaT299
    #print len(thetaT2prof)
    #print len(T2profAvg)


#    print pS
#    print pSArray
    np.savetxt(foutv2,np.c_[	dSArray,y,UeArray,deltasArray,thetaArray, \
		    		deltaT1WArray,deltaT2WArray,deltaT1Array,deltaT2Array, pSArray, \
				T1eArray, T1SArray,T2eArray,T2SArray,utauArray,\
				UprofAvg,VprofAvg, WprofAvg,T1profAvg,T2profAvg,\
				thetaT1prof, thetaT2prof, T1rmsprofAvg,T2rmsprofAvg,TUXprofAvg,\
				TUYprofAvg,TUZprofAvg,UVprofAvg,UWprofAvg,VWprofAvg],\
                            fmt="%15.12f "*30)
			    
			    #%15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f \
                            #     %15.12f %15.12f %15.12f %15.12f %15.12f\
                            #     %15.12f %15.12f %15.12f %15.12f %15.12f \
                            #     %15.12f %15.12f %15.12f  \
                            #     %15.12f %15.12f %15.12f %15.12f %15.12f")




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
    TuxEdge = TUXprofAvg[U99Index]
    TuyEdge = TUYprofAvg[U99Index]
    TuzEdge = TUZprofAvg[U99Index]

    UVEdge = UVprofAvg[U99Index]
    UWEdge = UWprofAvg[U99Index]
    VWEdge = VWprofAvg[U99Index]

    TuxMax = max(TUXprofAvg)
    TuyMax = max(TUYprofAvg)
    TuzMax = max(TUZprofAvg)

    UVMax = min(UVprofAvg)
    UWMax = max(UWprofAvg)
    VWMax = max(VWprofAvg)


    TuxMaxIndex = TUXprofAvg.index(TuxMax)
    TuyMaxIndex = TUYprofAvg.index(TuyMax)
    TuzMaxIndex = TUZprofAvg.index(TuzMax)

    UVMaxIndex = UVprofAvg.index(UVMax)
    UWMaxIndex = UWprofAvg.index(UWMax)
    VWMaxIndex = VWprofAvg.index(VWMax)


    yTuxMax = y[TuxMaxIndex]
    yTuyMax = y[TuyMaxIndex]
    yTuzMax = y[TuzMaxIndex]
    yUVMax = y[UVMaxIndex]
    yUWMax = y[UWMaxIndex]
    yVWMax = y[VWMaxIndex]

    UVMaxNon = UVMax / Ue**2
    #thermal
    T1rmsMax = max(T1rmsprofAvg)
    T2rmsMax = max(T2rmsprofAvg)

    T1rmsMaxIndex = T1rmsprofAvg.index(T1rmsMax)
    T2rmsMaxIndex = T2rmsprofAvg.index(T2rmsMax)
    if T1rmsMax ==0:
        yT1rmsMax =0.0
    else:
        yT1rmsMax = y[T1rmsMaxIndex]

    if T2rmsMax ==0:
        yT2rmsMax =0.0
    else:
        yT2rmsMax = y[T2rmsMaxIndex]

    #print "UVIndex = " + str(UVMaxIndex)
    #print "TuxIndex = " + str(TuxMaxIndex)
    #print "UVmax = " + str(UVMax)
    #print "TuxMax = " + str(TuxMax)
    #print "y @UVMAx = " +str(yUVMax)
    #print "y @TUxMax = " +str(yTuxMax)
    
    #print [x for x in UprofAvg if x > 0.99*Ue]
    #print [i for i in range(len(UprofAvg)) if UprofAvg[i] > 0.99*Ue]
    #print y[0:UeIndex+1]


    foutv11.write( "%15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f \
                   %15.12f %15.12f %15.12f %15.12f %15.12f  \
                    %15.12f %15.12f \
                    %15.12f %15.12f \
                    %15.12f %15.12f %15.12f %15.12f  \
                    %15.12f %15.12f  %15.12f %15.12f\
                   %15.12f %15.12f %15.12f %15.12f  \n" % \
                    (dS, Ue, pS, St1, St2, T1S, T2S, T1e, T2e, delta, \
                     delta95, delta99, delta999, deltas, theta,\
                      H, G, \
                     deltaT1, deltaT195, deltaT199, deltaT1999, \
                     deltaT2, deltaT295, deltaT299, deltaT2999, \
                     deltaT1W, deltaT2W, cf, utau ))

    str1="%15.12f "*29 + " \n"
    print str1
    foutv12.write( ("%15.12f "*29 +" \n") % \
                    (dS, Ue, T1S, T2S, T1e,\
		     T2e, TuxMax, TuyMax, TuzMax, UVMax,\
		     UWMax, VWMax, UVMaxNon, TuxEdge, TuyEdge,\
		     TuzEdge, UVEdge, UWEdge, VWEdge, T1rmsMax,\
		     T2rmsMax,yTuxMax, yTuyMax, yTuzMax, yUVMax,\
		     yUWMax,yVWMax, yT1rmsMax, yT2rmsMax))

    index=index+1		
    #print "index=" + str(index)
    del y
    del UprofAvg
    del UprofAvgAbs
    del VprofAvg
    del WprofAvg
    del T1profAvg
    del T2profAvg
    del T1rmsprofAvg
    del T2rmsprofAvg
#    del UT1profAvg
#    del UT2profAvg
    del TUXprofAvg
    del TUYprofAvg
    del TUZprofAvg
    del UVprofAvg
    del UWprofAvg
    del VWprofAvg
    del dSArray
    del UeArray
    del utauArray

    del T1eArray
    del T2eArray

    del T1SArray
    del T2SArray
foutv11.close()
foutv12.close()
foutv2.close()
del foutv2



Delete(integrateVariables2)
del integrateVariables2
Delete(slice2)
del slice2
Delete(clip4)
del clip4
Delete(clip3)
del clip3
Delete(clip2)
del clip2


Delete(integrateVariables1)
del integrateVariables1
Delete(plotOverLine1)
del plotOverLine1
Delete(clip1)
del clip1
Delete(slice1)
del slice1

 
                                                                                                                                                                                                                                                                                                        
