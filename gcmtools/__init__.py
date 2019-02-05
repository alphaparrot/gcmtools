import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap
import netCDF4 as nc

class DimensionError(Exception):
    pass

class UnitError(Exception):
    pass

class DatafileError(Exception):
    pass

def show():
    plt.show()

def parse(file,variable,lat=None,lon=None):
    ncd=nc.Dataset(file,"r")
    variable = ncd.variables[variable][:]
    
    if lat:
        lt = ncd.variables[lat][:]
    elif "lat" in ncd.variables:
        lt = ncd.variables['lat'][:]
    elif "lt" in ncd.variables:
        lt = ncd.variables["lt"][:]
    elif "lats" in ncd.variables:
        lt = ncd.variables['lats'][:]
    elif "latitude" in ncd.variables:
        lt = ncd.variables['latitude'][:]
    elif "latitudes" in ncd.variables:
        lt = ncd.variables['latitudes'][:]
    else:
        raise DatafileError("Unknown datafile format; unsure how to extract latitude")
    
    if lon:
        ln = ncd.variables[lon][:]
    elif "lat" in ncd.variables:
        ln = ncd.variables['lon'][:]
    elif "lt" in ncd.variables:
        ln = ncd.variables["ln"][:]
    elif "lats" in ncd.variables:
        ln = ncd.variables['lons'][:]
    elif "latitude" in ncd.variables:
        ln = ncd.variables['longitude'][:]
    elif "latitudes" in ncd.variables:
        ln = ncd.variables['longitudes'][:]
    else:
        raise DatafileError("Unknown datafile format; unsure how to extract longitude")
    
    ncd.close()
    
    return ln,lt,variable

def make2d(variable,lat=None,lon=None,time=None,lev=None,ignoreNaNs=True):
    if ignoreNaNs:
        sumop = np.nansum
        meanop = np.nanmean
    else:
        sumop = np.sum
        meanop = np.mean
    if len(variable.shape)==2:
        return variable
    if time:
        try:
            variable=variable[time,:]
        except:
            raise UnitError("You have probably passed a float time to a variable with no "+
                            "information about what that means. You should pass an integer "+
                            "time index instead")
    elif time==None and len(variable.shape)>2:
        variable=meanop(variable,axis=0)
    elif time==0:
        variable=variable[time,:]
    if len(variable.shape)>2:
        if lev!=None:
            if type(lev)==int:
                variable=variable[lev,:]
            elif lev=="sum":
                variable=sumop(variable,axis=0)
            elif lev=="mean":
                variable=meanop(variable,axis=0)
            else:
                raise UnitError("Unknown level specification")
        elif lat!=None and lon==None:
            if type(lat)==int:
                variable=variable[:,lat,:]
            elif lat=="sum":
                variable=sumop(variable,axis=1)
            elif lat=="mean":
                variable=meanop(variable,axis=1)
            else:
                raise UnitError("Unknown latitude specification")
        elif lon!=None and lat==None:
            if type(lon)==int:
                variable=variable[:,:,lon]
            elif lon=="sum":
                variable=sumop(variable,axis=2)
            elif lon=="mean":
                variable=meanop(variable,axis=2)
            else:
                raise UnitError("Unknown longitude specification")
        else:
            raise DimensionError("Inappropriate or insufficient dimensional constraints")
    
    return variable
    

def spatialmath(variable,lat=None,lon=None,file=None,mean=True,time=None,
               ignoreNaNs=True,lev=None,radius=6.371e6):
    
    if ignoreNaNs:
        sumop = np.nansum
        meanop = np.nanmean
    else:
        sumop = np.sum
        meanop = np.mean
        
    if file:
        ln,lt,variable = parse(file,variable,lat=lat,lon=lon)
        
    else:
        if lat==None or lon==None:
            raise DimensionError("Need to provide latitude and longitude data")
        ln=lon
        lt=lat
    variable = make2d(variable,time=time,lev=lev,ignoreNaNs=ignoreNaNs)
    
    lt1 = np.zeros(len(lt)+1)
    lt1[0] = 90
    for n in range(0,len(lt)-1):
        lt1[n+1] = 0.5*(lt[n]+lt[n+1])
    lt1[-1] = -90
    dln = np.diff(ln)[0]
    ln1 = np.zeros(len(ln)+1)
    ln1[0] = -dln
    for n in range(0,len(ln)-1):
        ln1[n+1] = 0.5*(ln[n]+ln[n+1])
    ln1[-1] = 360.0-dln
    
    lt1*=np.pi/180.0
    ln1*=np.pi/180.0
    
    darea = np.zeros((len(lt),len(ln)))
    for jlat in range(0,len(lt)):
        for jlon in range(0,len(ln)):
            dln = ln1[jlon+1]-ln1[jlon]
            darea[jlat,jlon] = (np.sin(lt1[jlat])-np.sin(lt1[jlat+1]))*dln
    
    svar = variable*darea
    if mean:
        outvar = sumop(svar)/sumop(darea)
    else:
        outvar = sumop(svar) * radius**2
    
    return outvar

def wrap2d(var):
    newvar = np.zeros(np.array(var.shape)+np.array((0,1)))
    newvar[:,:-1] = var[:,:]
    newvar[:,-1] = var[:,0]
    return newvar

def pcolormesh(variable,lon=None,lat=None,projection=None,cmap="viridis",
         shading='Gouraud',norm=None,vmin=None,vmax=None,invertx=False,
         inverty=False,linthresh=None,linscale=None,gamma=None,bounds=None,
         symmetric=False,ncolors=256,**kwargs):
    
    if symmetric==True: #assumes zero is the midpoint
        if vmin!=None and vmax!=None:
            vmin=-max(abs(vmin,vmax))
            vmax= max(abs(vmin,vmax))
        elif vmin!=None:
            vmax=-vmin
        elif vmax!=None:
            vmin=-vmax
        else:
            vmax = np.nanmax(abs(variable))
            vmin = -vmax
            
    elif symmetric: # a midpoint is specified
        if vmin!=None and vmax!=None:
            vmin = symmetric - max(abs(vmin,vmax))
            vmax = 2*symmetric - vmin
        elif vmin!=None:
            vmax = 2*symmetric - vmin
        elif vmax!=None:
            vmin = 2*symmetric - vmax
        else:
            vmax = symmetric + np.nanmax(abs(variable-symmetric))
            vmin = 2*symmetric - vmax
    
    if norm=="Log":
        normalization=colors.LogNorm(vmin=vmin,vmax=vmax)
    elif norm=="SymLog":
        normalization=colors.SymLog(vmin=vmin,vmax=vmax,linthresh=linthresh,linscale=linscale)
    elif norm=="PowerLog":
        normalization=colors.PowerLog(gamma,vmin=vmin,vmax=vmax)
    elif norm=="Bounds":
        normalization=colors.BoundaryNorm(bounds=bounds,ncolors=ncolors)
    else:
        normalization=colors.Normalize(vmin=vmin,vmax=vmax)
    
    if type(lon)==type(None) or type(lat)==type(None):
        im = plt.pcolormesh(variable,norm=normalization,shading=shading,cmap=cmap)
        if inverty:
            plt.gca().invert_yaxis()
        if invertx:
            plt.gca().invert_xaxis()
        return im
    
    if projection:
        
        if len(lon.shape)==1:
            lon,lat = np.meshgrid(lon,lat)
            lon = wrap2d(lon)
            lon[:,-1] = lon[:,0]+360.0
            lat = wrap2d(lat)
        variable=wrap2d(variable)
        m=Basemap(projection=projection,**kwargs)
        im=m.pcolormesh(lon,lat,variable,cmap=cmap,shading=shading,norm=normalization,latlon=True)
        return m,im
    
    im=plt.pcolormesh(lon,lat,variable,cmap=cmap,shading=shading,norm=normalization,**kwargs)
    if inverty:
        plt.gca().invert_yaxis()
    if invertx:
        plt.gca().invert_xaxis()
    return im
    