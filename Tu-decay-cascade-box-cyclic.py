import numpy as np
import os
import sys
#sys.modules[__name__].__dict__.clear()
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
print(os.getcwd())
#path1=np.genfromtxt("/home/03624/ykanani/Repositories/paraview-scripts/localOutputFolder",dtype='str')
#print(path1)
##############################
############# 1.55/100000/0.711*38000/(15.8*(T1Mean-20))

## find source
sourceName = raw_input("Please enter source name: ")
caseName = sourceName # raw_input("Please enter case name: ")
print sourceName
source = FindSource(sourceName + '.foam')
print source



n= 100
xmin=-0.252 
xmax=0.1
Dz=input("Please enter the span")
Dy=0.3848484754562378
ratioz=1 #ratio of the boxDz tpo Dz
ratioy=1 
nu = 1.55/100000.0
U = 4
area = Dz*ratioz*Dy*ratioy
print('area=',area)
#############
###############################


# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# create a new 'Slice'
slice1 = Slice(Input=source)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [0, 0, 0]

# create a new 'Integrate Variables'
integrateVariables1 = IntegrateVariables(Input=slice1)


location =input("Please enter the location that you are running this script, 0:Marvericks,1:Office PC,2:BW :")
if location==0:
	path1="/home/03624/ykanani/tempparaview/"
elif location==1:
	path1='E:/PostProcess/pvout/'
elif location==2:
	path1='/u/sciteam/kanani/tempparaview/'

path2=".plt"


path = path1 + caseName + "_Tu_" +str(xmin) + "to" +str(xmax) + "_n" + str(n)  + path2
fout = open(path,'a')
#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
Tux=[]
Tuy=[]
Tuz=[]
Tu=[]
xx=[]

for i in range(0,n):
    print i
    if (n==1):
        x=xmin
    else:
        x = (xmax-xmin)/(n-1)*i+xmin
    
    
    print x
    slice1.SliceType.Origin = [x, 0, 0]
    disp = Show(integrateVariables1, renderView1)
    integrateVariables1.UpdatePipeline()
    

    pdi = servermanager.Fetch(integrateVariables1)
    
    ##print pdi.GetCellData().GetArray("U").GetValue(1)
    ###print pdi.GetBlock(0)
    #print pdi.GetBlock(0).GetBlock(0)
    #print pdi.GetBlock(0).GetBlock(0).GetBlock(0).GetPointData()
    UPrime2Mean = pdi.GetCellData().GetArray("UPrime2Mean")  
 
    Uxx = UPrime2Mean.GetValue(0)
    Uyy = UPrime2Mean.GetValue(1)
    Uzz = UPrime2Mean.GetValue(2)
    
    Tux.append(np.sqrt(Uxx/area/U**2))
    Tuy.append(np.sqrt(Uyy/area/U**2)) 
    Tuz.append(np.sqrt(Uzz/area/U**2)) 
    
    Tu.append(np.sqrt((Uxx+Uyy+Uzz)/area/3/U**2))
    xx.append(x)
    
print Tux
print xx
    
data=np.c_[xx,Tux,Tuy,Tuz,Tu]


#np.savetxt(fout,St,fmt="%10.5f ")
np.savetxt(fout,data,header="VARIABLES=dx Tux Tuy Tuz Tu\nZONE T=\"" + caseName +"\"",comments='',fmt="%15.12f %15.12f %15.12f %15.12f %15.12f")
fout.close()
del Tux
del Tuy
del Tuz
del Tu
del xx
