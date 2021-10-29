import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import netCDF4 as nc

from matplotlib.pyplot import *

try:
    del globals()["pcolormesh"]
except:
    pass
try:
    del globals()["colors"]
except:
    pass

import matplotlib.colors as colors

def xcolorbar(mappable,fontsize=None,ticksize=None,nticks=None,**kwargs):
    cbar = plt.colorbar(mappable,**kwargs)
    if fontsize is not None and "label" in kwargs:
        cbar.set_label(label=kwargs["label"],fontsize=fontsize)
    if ticksize is not None:
        cbar.ax.tick_params(labelsize=ticksize)
        if nticks is not None:
            from matplotlib import ticker
            tick_locator = ticker.MaxNLocator(nbins=7)
            cbar.locator = tick_locator
            cbar.update_ticks()
    return cbar

class _xBasemap(Basemap):
    def __init__(self,**kwargs):
        super(_xBasemap,self).__init__(**kwargs)
    def xcolorbar(self,mappable,fontsize=None,ticksize=None,nticks=None,**kwargs):
        cbar = super(_xBasemap,self).colorbar(mappable,**kwargs)
        if fontsize is not None and "label" in kwargs:
            cbar.set_label(label=kwargs["label"],fontsize=fontsize)
        if ticksize is not None:
            cbar.ax.tick_params(labelsize=ticksize)
        if nticks is not None:
            from matplotlib import ticker
            tick_locator = ticker.MaxNLocator(nbins=7)
            cbar.locator = tick_locator
            cbar.update_ticks()
        return cbar

class _Dataset:
    def __init__(self,filename):
        if filename[-3:]==".nc":
            self.body = nc.Dataset(filename,"r")
            self.variables=self.body.variables
        elif filename[-4:]==".npy":
            self.body = np.load(filename)
            self.variables = self.body.item()
        else:
            raise DatafileError("Unknown dataset format")
        
    def close(self):
        try:
            self.body.close()
        except:
            pass

class DimensionError(Exception):
    pass

class UnitError(Exception):
    pass

class DatafileError(Exception):
    pass

def parse(file,variable,lat=None,lon=None):
    ncd=_Dataset(file)
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

def make2d(variable,lat=None,lon=None,time=None,lev=None,ignoreNaNs=True,radius=6.371e6,
           latitudes=None,longitudes=None):
    """Compress a variable in two dimensions by slicing or averaging.

    Parameters
    ----------
    variable : numpy.ndarray
        The variable to operate on
    lat,lon,lev : int, str, optional
        Either an index on which to slice, or either of "sum" or "mean", indicating what
        should be done along that axis.
    time : int, optional
        The time index on which to slice. If unspecified, a time average will be returned.
    ignoreNaNs : bool, optional
        If set, will use NaN-safe numpy operators.
    radius : float, optional
        Planet radius in meters (only used for summation)
    latitudes: numpy.ndarray, optional
        Latitude array--required if lat is "mean", or if either lat or lon is "sum"
    longitudes: numpy.ndarray, optional
        Longitude array--required if lon is "mean" or if either lat or lon is "sum"
        
    Returns
    -------
    numpy.ndarray
        A 2-D array
    """
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
                variable=latsum(variable,latitudes,dlon=longitudes[1]-longitudes[0],
                                radius=radius)
            elif lat=="mean":
                variable=latmean(variable,latitudes)
            else:
                raise UnitError("Unknown latitude specification")
        elif lon!=None and lat==None:
            if type(lon)==int:
                variable=variable[:,:,lon]
            elif lon=="sum":
                newvar = np.zeros(variable.shape[:-1])
                gradlat = np.gradient(np.sin(latitudes*np.pi/180.))
                for lt in range(len(latitudes)):
                    newvar[:,lt]=lonsum(variable,longitudes,dsinlat=gradlat[lt],radius=radius)
                variable = newvar
            elif lon=="mean":
                variable=lonmean(variable,longitudes)
            else:
                raise UnitError("Unknown longitude specification")
        else:
            raise DimensionError("Inappropriate or insufficient dimensional constraints")
    
    return variable
    

def spatialmath(variable,lat=None,lon=None,file=None,mean=True,time=None,
               ignoreNaNs=True,lev=None,radius=6.371e6):
    """Compute spatial means or sums of data

    Parameters
    ----------
    variable : str, numpy.ndarray
        The variable to operate on. Can either be a data array, or the name of a variable. If the latter, file must be specified.
    lat,lon : numpy.ndarray, optional
        Latitude and longitude arrays. If file is provided and lat and lon are not, they will be
        extracted from the file.
    file : str, optional
        Path to a NetCDF output file to open and extract data from.
    mean : bool, optional
        If True, compute a global mean. If False, compute a global sum.
    time : int, optional
        The time index on which to slice. If unspecified, a time average will be returned.
    ignoreNaNs : bool, optional
        If True, use NaN-safe numpy operators.
    lev : int, optional
        If set, slice a 3D spatial array at the specified level.
    radius : float, optional
        Radius of the planet in meters. Only used if mean=False.
        
    Returns
    -------
    float

    """
    
    if ignoreNaNs:
        sumop = np.nansum
        meanop = np.nanmean
    else:
        sumop = np.sum
        meanop = np.mean
        
    if file:
        ln,lt,variable = parse(file,variable,lat=lat,lon=lon)
        
    else:
        if type(lat)==type(None) or type(lon)==type(None):
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
            darea[jlat,jlon] = abs(np.sin(lt1[jlat])-np.sin(lt1[jlat+1]))*abs(dln)
    
    svar = variable*darea
    if mean:
        outvar = sumop(svar)/sumop(darea)
    else:
        outvar = sumop(svar) * radius**2
    
    return outvar

def latmean(variable,latitudes):
    """Compute meriodional mean (i.e. the variable that changes is latitude).
    
    Compute the area-weighted mean of a latitude array :math:`x`\ , such that:

    .. math::

        \\bar{x} = \\frac{\sum_{i=1}^N |\\sin(\\phi_{i-1/2})-\\sin(\\phi_{i+1/2})|x_i}{\sum_{i=1}^N |\\sin(\\phi_{i-1/2})-\\sin(\\phi_{i+1/2})|}
    
    Parameters
    ----------
    variable : numpy.ndarray
        Array to be averaged. Assumption is that if 2D, lat is the first dimension, if 3D, the second dimension, and if 4D. the 3rd dimension.
    latitudes : array-like
        Array or list of latitudes
        
    Returns
    -------
    scalar or numpy.ndarray
        Depending on the dimensionality of the input array, output may have 0, 1, or 2 dimensions.
    """
        
    lt1 = np.zeros(len(latitudes)+1)
    lt1[0] = 90
    for n in range(0,len(latitudes)-1):
        lt1[n+1] = 0.5*(latitudes[n]+latitudes[n+1])
    lt1[-1] = -90
    
    lt1*=np.pi/180.0
    
    darea = np.zeros(latitudes.shape)
    for jlat in range(0,len(latitudes)):
        darea[jlat] = np.sin(lt1[jlat])-np.sin(lt1[jlat+1])
    
    if len(variable.shape)==1:
        return np.nansum(variable*darea)/np.nansum(darea)
    elif len(variable.shape)==2:
        return np.nansum(variable*darea[:,np.newaxis],axis=0)/np.nansum(darea[:,np.newaxis]*np.ones(variable.shape),axis=0)
    elif len(variable.shape)==3:
        return np.nansum(variable*darea[np.newaxis,:,np.newaxis],axis=1)/np.nansum(darea[np.newaxis,:,np.newaxis]*np.ones(variable.shape),axis=1)
    elif len(variable.shape)==4:
        return np.nansum(variable*darea[np.newaxis,np.newaxis,:,np.newaxis],axis=2)/np.nansum(darea[np.newaxis,np.newaxis,:,np.newaxis]*np.ones(variable.shape),axis=2)
    else:
        raise DimensionError("Variable must have 4 or fewer dimensions. Latitude should be the second-from the right-most dimension if there are 2 or more dimensions.")
        
def latsum(variable,latitudes,dlon=360.0,radius=6.371e6):
    """Compute meriodional sum (i.e. the variable that changes is latitude).
    
    Compute the area-weighted sum of a latitude array :math:`x` given a longitude span :math:`\\Delta\\theta` and planet radius :math:`R`\ , such that:

    .. math::

        X = \sum_{i=1}^N |\\sin(\\phi_{i-1/2})-\\sin(\\phi_{i+1/2})|\\Delta\\theta R^2x_i
    
    Parameters
    ----------
    variable : numpy.ndarray
        Array to be summed. Assumption is that if 2D, lat is the first dimension, if 3D, the second dimension, and if 4D. the 3rd dimension.
    latitudes : array-like
        Array or list of latitudes
    dlon : float, optional
        Longitude span in degrees.
    radius : float, optional
        Planet radius in meters.
        
    Returns
    -------
    scalar or numpy.ndarray
        Depending on the dimensionality of the input array, output may have 0, 1, or 2 dimensions.
    """
        
    lt1 = np.zeros(len(latitudes)+1)
    lt1[0] = 90
    for n in range(0,len(latitudes)-1):
        lt1[n+1] = 0.5*(latitudes[n]+latitudes[n+1])
    lt1[-1] = -90
    
    lt1*=np.pi/180.0
    dlon *= np.pi/180.0
    
    darea = np.zeros(latitudes.shape)
    for jlat in range(0,len(latitudes)):
        darea[jlat] = abs(np.sin(lt1[jlat])-np.sin(lt1[jlat+1]))*abs(dlon)*radius**2
    
    if len(variable.shape)==1:
        return np.nansum(variable*darea)
    elif len(variable.shape)==2:
        return np.nansum(variable*darea[:,np.newaxis],axis=0)
    elif len(variable.shape)==3:
        return np.nansum(variable*darea[np.newaxis,:,np.newaxis],axis=1)
    elif len(variable.shape)==4:
        return np.nansum(variable*darea[np.newaxis,np.newaxis,:,np.newaxis],axis=2)
    else:
        raise DimensionError("Variable must have 4 or fewer dimensions. Latitude should be the second-from the right-most dimension if there are 2 or more dimensions.")


def lonmean(variable,longitudes):
    """Compute zonal mean (i.e. the variable that changes is longitude).
    
    Compute the area-weighted mean of a longitude array :math:`x`\ , such that:

    .. math::

        \\bar{x} = \\frac{\sum_{i=1}^N |\\theta_{i-1/2}-\\theta_{i+1/2}|x_i}{\sum_{i=1}^N |\\theta_{i-1/2}-\\theta_{i+1/2}|}
    
    Parameters
    ----------
    variable : numpy.ndarray
        Array to be summed. Assumption is that longitude is always the last dimension.
        
    Returns
    -------
    scalar or numpy.ndarray
        Depending on the dimensionality of the input array, output may be a scalar or have N-1 dimensions.
    """
    
    dlon = np.gradient(longitudes)
    sumlon = np.nansum(dlon)
    dlon = np.broadcast_to(dlon,variable.shape)
    
    return np.nansum(variable*dlon,axis=-1)/sumlon
        
def lonsum(variable,longitudes,dsinlat=2.0,radius=6.371e6):
    """Compute zonal sum (i.e. the variable that changes is longitude).
    
    Compute the area-weighted sum of a longitude array :math:`x` given a latitude span :math:`\\Delta\\sin\\phi` and planet radius :math:`R`\ , such that:

    .. math::

        X = \sum_{i=1}^N |\\theta_{i-1/2}-\\theta_{i+1/2}|\\Delta\\sin\\phi R^2x_i
    
    Parameters
    ----------
    variable : numpy.ndarray
        Array to be summed. Assumption is that longitude is always the last dimension.
    longitudes : array-like
        Array or list of longitudes
    dsinlat : float, optional
        The sine-latitude span for the longitude span considered. The default is 2, corresponding to -90 degrees to 90 degrees.
    radius : float, optional
        Planet radius in meters.
        
    Returns
    -------
    scalar or numpy.ndarray
        Depending on the dimensionality of the input array, output may have 0, 1, or 2 dimensions.
    """
        
    dlon = np.gradient(longitudes)*np.pi/180.0
    
    darea = np.zeros(longitudes.shape)
    for jlon in range(0,len(longitudes)):
        darea[jlon] = abs(dsinlat)*abs(dlon[jlon])*radius**2
    
    darea = np.broadcast_to(darea,variable.shape)
    
    return np.nansum(variable*darea,axis=-1)
    
def cspatialmath(variable,lat=None,lon=None,file=None,mean=True,time=None,
               ignoreNaNs=True,lev=None,radius=6.371e6,poles=False):
    """Compute spatial means or sums of data, but optionally don't go all the way to the poles.

    Sometimes, saying that the latitudes covered go all the way to :math:`\pm90^\circ` results in
    errors, and accurate accounting requires excluding the poles themselves. This function
    is identical to spatialmath, except that it provides that option.

    Parameters
    ----------
    variable : str, numpy.ndarray
        The variable to operate on. Can either be a data array, or the name of a variable. If the latter, file must be specified.
    lat,lon : numpy.ndarray, optional
        Latitude and longitude arrays. If file is provided and lat and lon are not, they will be
        extracted from the file.
    file : str, optional
        Path to a NetCDF output file to open and extract data from.
    mean : bool, optional
        If True, compute a global mean. If False, compute a global sum.
    time : int, optional
        The time index on which to slice. If unspecified, a time average will be returned.
    ignoreNaNs : bool, optional
        If True, use NaN-safe numpy operators.
    lev : int, optional
        If set, slice a 3D spatial array at the specified level.
    radius : float, optional
        Radius of the planet in meters. Only used if mean=False.
    poles : bool, optional
        If False (default), exclude the poles.
        
    Returns
    -------
    float

    """
    
    if ignoreNaNs:
        sumop = np.nansum
        meanop = np.nanmean
    else:
        sumop = np.sum
        meanop = np.mean
        
    if file:
        ln,lt,variable = parse(file,variable,lat=lat,lon=lon)
        
    else:
        if type(lat)==type(None) or type(lon)==type(None):
            raise DimensionError("Need to provide latitude and longitude data")
        ln=lon
        lt=lat
    variable = make2d(variable,time=time,lev=lev,ignoreNaNs=ignoreNaNs)
    
    lt1 = np.zeros(len(lt)+1)
    dlt1 = abs(np.diff(lt)[0])
    dlt2 = abs(np.diff(lt)[-1])
    lt1[0] = lt[0]+0.5*dlt1
    for n in range(0,len(lt)-1):
        lt1[n+1] = 0.5*(lt[n]+lt[n+1])
    lt1[-1] = lt[-1]-0.5*dlt2
    if poles:
        lt1[0] = 90.0
        lt1[-1] = -90.0
    
    dln = np.diff(ln)[0]
    ln1 = np.zeros(len(ln)+1)
    ln1[0] = ln[0]-dln*0.5
    for n in range(0,len(ln)-1):
        ln1[n+1] = 0.5*(ln[n]+ln[n+1])
    ln1[-1] = ln[-1]+dln*0.5
    
    lt1*=np.pi/180.0
    ln1*=np.pi/180.0
    
    darea = np.zeros((len(lt),len(ln)))
    for jlat in range(0,len(lt)):
        for jlon in range(0,len(ln)):
            dln = ln1[jlon+1]-ln1[jlon]
            darea[jlat,jlon] = abs(np.sin(lt1[jlat])-np.sin(lt1[jlat+1]))*abs(dln)
    
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
         inverty=False,linthresh=1.0e-3,linscale=1.0,gamma=1.0,bounds=None,
         symmetric=False,ncolors=256,**kwargs):
    
    if symmetric==True: #assumes zero is the midpoint
        if vmin and vmax:
            vmin=-max(abs(vmin),abs(vmax))
            vmax= max(abs(vmin),abs(vmax))
        elif vmin:
            vmax=-vmin
        elif vmax:
            vmin=-vmax
        else:
            vmax = np.nanmax(abs(variable))
            vmin = -vmax
            
    elif type(symmetric)!=type(None) and type(symmetric)!=type(False): # a midpoint is specified
        if vmin and vmax:
            vmin = symmetric - max(abs(vmin-symmetric),abs(vmax-symmetric))
            vmax = 2*symmetric - vmin
        elif vmin:
            vmax = 2*symmetric - vmin
        elif vmax:
            vmin = 2*symmetric - vmax
        else:
            vmax = symmetric + np.nanmax(abs(variable-symmetric))
            vmin = 2*symmetric - vmax
    
    else:
        if not vmin:
            vmin = np.nanmin(variable)
        if not vmax:
            vmax = np.nanmax(variable)
    
    if norm=="Log":
        normalization=colors.LogNorm(vmin=vmin,vmax=vmax)
    elif norm=="SymLog":
        normalization=colors.SymLogNorm(vmin=vmin,vmax=vmax,linthresh=linthresh,linscale=linscale)
    elif norm=="PowerLog":
        normalization=colors.PowerNorm(gamma,vmin=vmin,vmax=vmax)
    elif norm=="Bounds":
        if type(bounds)==type(None):
            bounds = np.linspace(vmin,vmax,num=ncolors+1)
        normalization=colors.BoundaryNorm(bounds,ncolors)
    else:
        normalization=colors.Normalize(vmin=vmin,vmax=vmax)
    
    if len(variable.shape)>2:
        variable=make2d(variable)
    
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
        m=_xBasemap(projection=projection,**kwargs)
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

def hadley(file,time=None,contours=None,ylog=False,**kwargs):
    ln,lt,levs=parse(file,"lev")
    plevs = levs*spatialmath("ps",file=file,time=time)
    ln,lt,va = parse(file,"va")
    va = make2d(va,lon="mean",time=time)
    strf = _stream(va,lt,plevs)
    ln,lt,ua = parse(file,"ua")
    uavg=make2d(ua,lon="mean",time=time)
    umin = 10*(int(np.nanmin(uavg))/10+1)
    umax = 10*(int(np.nanmax(uavg))/10)
    clvs = np.linspace(umin,umax,num=(umax-umin)/10+1)
    if len(clvs) < 10:
        umin = 5*(int(np.nanmin(uavg))/5+1)
        umax = 5*(int(np.nanmax(uavg))/5)
        clvs = np.linspace(umin,umax,num=(umax-umin)/5+1)
    im=pcolormesh(strf,x=lt,y=plevs,cmap='RdBu_r',symmetric=True,invertx=True,inverty=True,**kwargs)
    plt.colorbar(im,label="Streamfunction [kg/s]")
    plt.xlabel("Latitude [$^\circ$N]")
    plt.ylabel("Pressure [hPa]")
    if ylog:
        plt.yscale('log')
    if contours:
        cs=plt.contour(lt,plevs,uavg,np.linspace(-100,100,num=41),colors='gray',linestyles='-')
        plt.clabel(cs,clvs,fmt='%1d',fontsize=8)
        return im,cs
    else:
        return im

    