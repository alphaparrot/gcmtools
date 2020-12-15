# gcmtools

By Adiv Paradise

_DEPENDENCIES: matplotlib, numpy, basemap, python-netcdf4_
_OPTIONAL (for interactive plots): nodejs, jupyter-widgets, jupyter-matplotlib_

This is a collection of custom functions, wrappers, and other tools to make analyzing 
output from a gcm just slightly easier. Mostly it functions as a wrapper for matplotlib
and basemap, but it includes useful things like a function for area-weighted math and a
function for computing the streamfunction and plotting the Hadley cells.

## Installation

You can either download and build/use from this repository, or you can use pip:

`pip install gcmtools`


## Usage

### parse(filename,variable,**kwargs)
 
Returns the data contained in variable, along with the latitude and longitude arrays.
 
    * filename
    
        Specifies the name of the netCDF4 file to open
        
    * variable
    
        Specifies the variable name to use
        
    * lat (optional)
    
        What name to use for the latitude array when parsing the file
        
    * lon (optional)
    
        What name to use for the longitude array when parsing the file
    
    
### make2d(variable,**kwargs)

Returns a 2D slice of the given variable
    
    * variable
    
        The data array to slice
        
    * ignoreNaNs (optional)
    
        Ignore NaNs when doing arithmetic operations (default=True)
        
    * lat (optional)
    
        type(int): slice the array at this latitude
        "sum": Take the meridional sum
        "mean": Take the meridional mean
        
    * lon (optional)
    
        type(int): slice the array at this longitude
        "sum" Take the zonal sum
        "mean": Take the zonal mean
        
    * lev (optional)
    
        type(int): slice the array at this vertical level
        "sum" Take the column sum
        "mean": Take the column mean
        
    * time (optional)
    
        type(int): Take the snapshot of the array at this time
        None (default): Take the time-average of the data
            

### spatialmath(variable,**kwargs)
    
Returns the area-weighted average or sum of the given variable
    
    * variable
    
        Either the name of the variable to use, or the data array itself.
        If the file keyword is used, this should be the variable name. If
        not, then the lat and lon arrays must be provided.
        
    * file (optional)
    
        The name of the file from which to extract the data
        
    * lat (optional)
    
        The latitude array to use (ignored if file keyword is used)
        
    * lon (optional)
    
        The longitude array to use (ignored if file keyword is used)
        
    * lev (optional)
    
        The level slice to use (see make2d() keyword options)
        
    * time (optional)
    
        type(int): Use the snapshot of the variable at this time
        None (default): Use the time-average of the variable
        
    * mean (optional)
    
        If True (default), the global mean will be calculated. If False,
        only the global sum will be returned.
        
    * radius (optional)
    
        The physical radius of the sphere with which to scale the sum 
        (if not computing the mean)

            
### wrap2d(variable)

Add a longitude column to a 2D lat-lon array, and fill it with the first column
    

### pcolormesh(variable,**kwargs)
    
Create and return a pcolormesh object showing variable. **kwargs can include all
normal pcolormesh keyword arguments, and if the 'projection' keyword argument is
specified, **kwargs can also contain any Basemap arguments.
    
gcmtools-specific arguments:
    
    * invertx
    
        Invert the x-axis. This is analogous to plt.gca().invert_xaxis()
        
    * inverty
    
        Invert the y-axis. This is analogous to plt.gca().invert_yaxis()
        
    * symmetric
    
        If True, compute a colormap normalization which is symmetric about zero.
        If not None and equal to a number, compute a colormap normalization 
        symmetric about that number. Useful for divergent colormaps.
    
Example:

`pcolormesh(temperature,x=lons,y=lats,projection='moll',lon_0=0,cmap='RdBu_r',symmetric=273.15)`

### hadley(filename,**kwargs)

Compute the streamfunction, and plot the zonal mean as a function of latitude and pressure. Optionally overplot zonal wind contours.

    * filename

       File from which to extract the streamfunction.
       
    * contours (optional)

       If True, compute the mean zonal wind and overplot it as a series of labeled contours.
       
    * ylog (optional)
    
       If True, use a logarithmic scale on the y-axis (corresponding to being linear in altitude).
       
       