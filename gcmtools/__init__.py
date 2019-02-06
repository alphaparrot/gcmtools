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

def pcolormesh(variable,x=None,y=None,projection=None,cmap="viridis",
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
    
    if type(x)==type(None) or type(y)==type(None):
        im = plt.pcolormesh(variable,norm=normalization,shading=shading,cmap=cmap)
        if inverty:
            plt.gca().invert_yaxis()
        if invertx:
            plt.gca().invert_xaxis()
        return im
    
    if projection:
        
        if len(x.shape)==1:
            x,y = np.meshgrid(x,y)
            x = wrap2d(x)
            x[:,-1] = x[:,0]+360.0
            y = wrap2d(y)
        variable=wrap2d(variable)
        m=Basemap(projection=projection,**kwargs)
        im=m.pcolormesh(x,y,variable,cmap=cmap,shading=shading,norm=normalization,latlon=True)
        return m,im
    
    im=plt.pcolormesh(x,y,variable,cmap=cmap,shading=shading,norm=normalization,**kwargs)
    if inverty:
        plt.gca().invert_yaxis()
    if invertx:
        plt.gca().invert_xaxis()
    return im
    
def _stream(va,lt,plevs):
    '''assumed that va has a shape of ([plevs],[lt])'''
    strf = np.zeros(va.shape)
    pref = 2*np.pi*6.371e6*np.cos(lt*np.pi/180.0)/9.81
    ps = np.array([0.0,]+list(plevs))
    vas = np.zeros(np.array(va.shape)+np.array((1,0)))
    vas[1:,:] = va[:,:]
    for k in range(0,len(plevs)):
        strf[k,:] = pref[:]*np.trapz(vas[0:k+1,:],x=ps[0:k+1],axis=0)
    return strf

def hadley(file,time=None,contours=None,ylog=False):
    ln,lt,levs=parse(file,"lev")
    plevs = levs*spatialmath("ps",file=file,time=time)
    ln,lt,va = parse(file,"va")
    va = make2d(va,lon="mean",time=time)
    strf = _stream(va,lt,plevs)
    ln,lt,ua = parse(file,"ua")
    uavg=make2d(ua,lon="mean",time=time)
    umin = 10*(int(uavg.min())/10+1)
    umax = 10*(int(uavg.max())/10)
    clvs = np.linspace(umin,umax,num=(umax-umin)/10+1)
    if len(clvs) < 10:
        umin = 5*(int(uavg.min())/5+1)
        umax = 5*(int(uavg.max())/5)
        clvs = np.linspace(umin,umax,num=(umax-umin)/5+1)
    im=pcolormesh(strf,x=lt,y=plevs,cmap='RdBu_r',symmetric=True,invertx=True,inverty=True)
    if contours:
        cs=plt.contour(lt,plevs,uavg,np.linspace(-100,100,num=41),colors='gray',linestyles='-')
        plt.clabel(cs,clvs,fmt='%1d',fontsize=8)
    plt.colorbar(im,label="Streamfunction [kg/s]")
    plt.xlabel("Latitude [$^\circ$N]")
    plt.ylabel("Pressure [hPa]")
    if ylog:
        plt.yscale('log')
    return im,cs
        
def savefig(filename,**kwargs):
    plt.savefig(filename,**kwargs)
    
def xlabel(label,**kwargs):
    plt.xlabel(label,**kwargs)
    
def ylabel(label,**kwargs):
    plt.ylabel(label,**kwargs)
    
def xscale(scaling,**kwargs):
    plt.xscale(scaling,**kwargs)
    
def yscale(scaling,**kwargs):
    plt.yscale(scaling,**kwargs)
    
def title(title,**kwargs):
    plt.title(title,**kwargs)
    
def subplots(*args,**kwargs):
    f,a=plt.subplots(*args,**kwargs)
    return f,a
    