"""Author: Ruth Green
   Calculate gradients for scalars, vectors, products of vectors
   Assumes dimensions lat and lon in degrees, pfull in hPa"""

import numpy as np
import xarray as xr
#from finite_difference import cfd

############
# mek236's code from execlim/sharecode/atmosphere/finite_difference.py
############

def cfd(F_in, x_in, axis=0, inc_bdy='both', cyclic=False):
    """
    Vectorized central finite difference for N-D array F_in by x_in.
    Parameters
    ----------
    F_in : array_like, N-dimensions
            Array to be differentiated
    x_in : array_like, either 1-D or F_in.shape == x_in.shape
            Coordinate array
    axis : integer
            Dimension where F_in.shape[axis] == x_in.shape[axis]
    inc_bdy : String, optional
            Include 0 'low' or N 'high', or both boundary(ies) in calculations
    cyclic : Boolean, optional
            Data is cyclic on <axis> if true
    Returns
    -------
    Centred
    """
    F = np.swapaxes(F_in, axis, 0)
    if len(x_in.shape) == len(F_in.shape):
        x = np.swapaxes(x_in, axis, 0)
    else:
        x = x_in.copy()
        for i in range(len(F_in.shape) - len(x_in.shape)):
            x = x[:, None]

    if inc_bdy in ['low', 'Low', 'L', 'l']:
        inc_l = True
        inc_h = False
    elif inc_bdy in ['high', 'High', 'H', 'h']:
        inc_l = False
        inc_h = True
    else:
        inc_l = True
        inc_h = True

    dFdx = np.zeros(F.shape)
    dFdx[1:-1, ...] = (F[2:, ...] - F[:-2, ...])/(x[2:, ...] - x[:-2, ...])
    if cyclic:
        dFdx[0, ...] = (F[-1, ...] - F[1, ...])/(x[-1, ...] - x[1, ...])
        dFdx[-1, ...] = (F[-2, ...] - F[0, ...])/(x[-2, ...] - x[0, ...])
    else:
        if inc_l:
            dFdx[0, ...] = (F[0, ...] - F[1, ...])/(x[0, ...] - x[1, ...])
        if inc_h:
            dFdx[-1, ...] = (F[-2, ...] - F[-1, ...])/(x[-2, ...] - x[-1, ...])

    return np.swapaxes(dFdx, axis, 0)

##################
##################


def ddx(field, a = 6376.0e3):
    """Calculate d/dx of a given DataArray. DataArray must include lat and lon dimensions"""
    
    try:
        field.coords['lon']
    except:
        raise NameError('Coord lon not found')
    try:
        field.coords['lat']
    except:
        raise NameError('Coord lat not found')
    
    coslat = np.cos(field.lat * np.pi/180)
    field_dx = cfd( field.values, field.lon.values*np.pi/180, field.get_axis_num('lon') )   
    field_dx = xr.DataArray( field_dx, dims = field.dims, coords = field.coords )
    field_dx = field_dx/coslat/a
    
    return field_dx


def ddy(field, vector = True, uv = False, a = 6376.0e3):
    """Calculate d/dy of a given DataArray. DataArray must include lat dimension.
        kwargs: vector - specify if input field is vector or scalar
                prod   - if a vector, is the field uv?"""
    
    try:
        field.coords['lat']
    except:
        raise NameError('Coord lat not found')
        
    coslat = np.cos(field.lat * np.pi/180)
    
    if vector and uv:
        cosfac = coslat**2
    elif vector:
        cosfac = coslat
    else:
        cosfac = 1.
    
    field_dy = cfd( (field*cosfac).values, field.lat.values*np.pi/180, field.get_axis_num('lat') )   
    field_dy = xr.DataArray( field_dy, dims = field.dims, coords = field.coords )
    field_dy = field_dy/cosfac/a
    
    return field_dy
    
    
def ddp(field):
    """Calculate d/dp of a given DataArray. DataArray must include pfull dimension"""
    
    try:
        field.coords['pfull']
    except:
        raise NameError('Coord pfull not found')
    
    field_dp = cfd( field.values, field.pfull.values*100., field.get_axis_num('pfull') )   
    field_dp = xr.DataArray( field_dp, dims = field.dims, coords = field.coords )
    
    return field_dp


def ddt(field, timedir = 'xofyear', secperunit = 5.*86400.):
    """Calculate d/dt in unit/s of a given DataArray. DataArray must include a time dimension
       Define seconds per unit time using secperunit. Default calc is for pentads"""
    
    try:
        field.coords[timedir]
    except:
        raise NameError('Coord ' + timedir + ' not found')
    
    field_dt = cfd( field.values, field.coords[timedir].values*secperunit, field.get_axis_num(timedir) )   
    field_dt = xr.DataArray( field_dt, dims = field.dims, coords = field.coords )
    
    return field_dt



if __name__ == '__main__':
    """Examples/sanity check"""
    import matplotlib.pyplot as plt
    from data_handling import time_means
    
    data = time_means('full_qflux', [121,481], filename='plev_pentad', timeav='pentad')
    
    dudx = ddx(data.ucomp[40,17,:,:])
    dudy = ddy(data.ucomp[40,17,:,:])
    dudp = ddp(data.ucomp[40,:,:,:].mean('lon'))
    dudt = ddt(data.ucomp[:,17,32,0])
    
    data.ucomp[40,17,:,:].plot.contourf(levels=np.arange(-50.,51.,10.))
    plt.figure(2)
    dudx.plot.contourf(levels=np.arange(-0.000016,0.0000161,0.000002))
    plt.figure(3)
    dudy.plot.contourf(levels=np.arange(-0.000016,0.0000161,0.000002))
    
    plt.figure(4)
    data.ucomp[40,:,:,:].mean('lon').plot.contourf(levels=np.arange(-50.,51.,10.), yincrease=False)
    plt.figure(5)
    dudp.plot.contourf(levels=np.arange(-0.0045,0.0046,0.0005), yincrease=False)
    
    plt.figure(6)
    data.ucomp[:,17,32,0].plot()
    plt.figure(7)
    dudt.plot()
    
    plt.show()