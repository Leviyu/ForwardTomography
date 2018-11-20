#!/usr/bin/python2.7
from numpy import *
from evtk.hl import gridToVTK 
from pylab import *
from sys import argv
#from mpl_toolkits.basemap import Basemap, cm


# INPUT:
#   model_file: input tomography xyz file
#   out_file:   output tomography vtk fike
#   nx, ny, nz: number of dep/lat/lon dimensions
#



script, model_file, out_file, nx,ny,nz= argv

nx = int(nx)
ny = int(ny)
nz = int(nz)
#//print model_file,nz

print script, model_file,out_file
print "nx ny nz", nx, ny, nz

close('all')
datafile=model_file

data=loadtxt(datafile)
x = data[:,0]
y = data[:,1]
z = data[:,2]
vs= data[:,3]
dep= data[:,4]
lat= data[:,5]
lon= data[:,6]


##nx = 30
##ny = 181
##nz = 361
print "Total line", nx*ny*nz

xx=x.reshape(nx,ny,nz)
yy=y.reshape(nx,ny,nz)
zz=z.reshape(nx,ny,nz)
dvsdvs=vs.reshape((nx,ny,nz))
depdep=dep.reshape((nx,ny,nz))
latlat=lat.reshape((nx,ny,nz))
lonlon=lon.reshape((nx,ny,nz))


##x = zeros((ndepth,nlon,nlat))
##y = zeros((ndepth,nlon,nlat))
##z = zeros((ndepth,nlon,nlat))
##XX = zeros((ndepth,nlon,nlat))
##YY = zeros((ndepth,nlon,nlat))
##ZZ = zeros((ndepth,nlon,nlat))
##new_dvs = zeros((ndepth,nlon,nlat))
XX = zeros((nx,ny,nz))
YY = zeros((nx,ny,nz))
ZZ = zeros((nx,ny,nz))
DVS = zeros((nx,ny,nz))
DEPDEP = zeros((nx,ny,nz))
LATLAT = zeros((nx,ny,nz))
LONLON = zeros((nx,ny,nz))

for i in range(nx): 
    for j in range(ny):
        for k in range(nz): 
            DVS[i,j,k] = dvsdvs[i,j,k]
            XX[i,j,k] = xx[i,j,k]
            YY[i,j,k] = yy[i,j,k]
            ZZ[i,j,k] = zz[i,j,k]
            DEPDEP[i,j,k] = depdep[i,j,k]
            LATLAT[i,j,k] = latlat[i,j,k]
            LONLON[i,j,k] = lonlon[i,j,k]


            ##print XX[i,j,k],YY[i,j,k], ZZ[i,j,k], DVS[i,j,k]
            ##YY[i,j,k] = (6371-DEPTH[i])*cos(LAT[k])*sin(LON[j])
            ##ZZ[i,j,k] = (6371-DEPTH[i])*sin(LAT[k])
            ##print XX[i,j,k]
            ##print x[i,j,k],y[i,j,k],z[i,j,k], new_dvs[i,j,k]

            ##new_dvs[i,j,k] = DVS[k,j,i]
            ##print x[i,j,k],y[i,j,k],z[i,j,k], vsarray[i,j,k]
            ##x[i,j,k] = DEPTH[i]
            ##y[i,j,k] = LON[j] 
            ##z[i,j,k] = LAT[k] 
##print "xx shape",xx.shape
##print "XX shape",XX.shape

##print(vsarray.shape)
##print "XX shape", XX.shape
##print "vs shape", vsarray.shape

## groofing around to make it work as well with python 3
##NumPy_data_shape = vsarray.shape
##VTK_data = numpy_support.numpy_to_vtk(num_array=vsarray.ravel(), deep=True, array_type=vtk.VTK_FLOAT)

# writing vts file, works with Python 2.X

##gridToVTK("./gypsum_test", XX, YY,ZZ, pointData = {"vs" : vsarray, "depth" : adepth, "lon" : alon, "lat": alat})
gridToVTK(out_file,XX,YY,ZZ, pointData = {"dvs" : DVS ,"depth" : DEPDEP ,"lat" : LATLAT, "Lon" : LONLON})










