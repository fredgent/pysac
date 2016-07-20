# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 17:45:48 2014

@author: sm1fg
"""

import numpy as np
import os
import astropy.units as u
import astropy
l_astropy_0=True
if astropy.__version__[0]=='1':
        l_astropy_0=False

hmi_model = {'photo_scale': 0.6*u.Mm,       #scale height for photosphere
             'chrom_scale': 0.1*u.Mm,      #scale height for chromosphere
             'corona_scale': 2.5e3*u.Mm,      #scale height for the corona
             'coratio': 0*u.one,  #footpoint portion scaling as corona 
             'model': 'hmi_model',
             'phratio': 0*u.one,  #footpoint portion scaling as photosphere
             'pixel': 0.36562475*u.Mm,      #(HMI pixel)
             'radial_scale': 0.10979002*u.Mm, #=> FWHM = half pixel
             'nftubes': 1,
             'B_corona': 0.*u.T,
             'pBplus': 1e-3*u.T}
hmi_model['chratio'] = 1*u.one - hmi_model['coratio'] - hmi_model['phratio']

sunspot = {'photo_scale': 0.60*u.Mm,
             'chrom_scale': 0.1*u.Mm,
             'corona_scale': 2.5e3*u.Mm,         #scale height for the corona
             'coratio': 0.03*u.one,
             'model': 'sunspot',
             'phratio': 0.0*u.one,
             'pixel': 0.36562475*u.Mm,  #(HMI pixel)
             'radial_scale': 3.6096*u.Mm,
             #'radial_scale': 0.14*u.Mm,
             'nftubes': 1,
             #'B_corona': 4.85e-4*u.T,
             'B_corona': 5.5e-4*u.T,
             'pBplus': 12.0e-4*u.T}
sunspot['chratio'] = 1*u.one - sunspot['coratio'] - sunspot['phratio']
#if 1D or 2D set unused dimensions to 0, and unrequired xyz limits to 1.
sunspot['Nxyz'] = [512,512,256] # 3D grid
#sunspot['Nxyz'] = [128,128,64] # 3D grid
sunspot['xyz']  = [-25*u.Mm,25*u.Mm,-25*u.Mm,25*u.Mm,2.5*u.Mm,22.5*u.Mm] #grid size

mfe_setup = {'photo_scale': 0.60*u.Mm,
             'chrom_scale': 0.4*u.Mm,
             'corona_scale': 0.25*u.Mm,  #scale height for the corona
             'coratio': 0.0*u.one,
             'model': 'mfe_setup',
             'phratio': 0.0*u.one,
             'pixel': 0.36562475*u.Mm,  #(HMI pixel)
             'radial_scale': 0.03938*u.Mm,
             #'radial_scale': 0.14*u.Mm,
             'nftubes': 1,
             #'B_corona': 4.85e-4*u.T,
             'B_corona': 5.5e-4*u.T,
             'pBplus': 12.0e-4*u.T}
mfe_setup['chratio'] = 1*u.one - mfe_setup['coratio'] - mfe_setup['phratio']
#if 1D or 2D set unused dimensions to 0, and unrequired xyz limits to 1.
mfe_setup['Nxyz'] = [128,128,128] # 3D grid
mfe_setup['Nxyz'] = [129,129,128] # 3D grid
mfe_setup['xyz']  = [-1*u.Mm,1*u.Mm,-1*u.Mm,1*u.Mm,3.6641221e-2*u.Mm,1.5877863*u.Mm] #grid size

spruit = {'photo_scale': 1.5*u.Mm,
          'chrom_scale': 0.5*u.Mm,
          'corona_scale': 100*u.Mm,      #scale height for the corona
          'coratio': 0.0*u.one,
          'model': 'spruit',
          'phratio': 0.0*u.one,
          'pixel': 0.1*u.Mm,              
          'radial_scale': 0.075*u.Mm,
          'nftubes': 1,
          'p0': 117200.0 * u.dyne/u.cm**2,
          'B_corona': 0.*u.T,
          'pBplus': 4.250e-4*u.T}
spruit['chratio'] = 1*u.one - spruit['coratio'] - spruit['phratio']
spruit['Nxyz'] = [128,128,256] # 3D grid
spruit['xyz']  = [-1.27*u.Mm,1.27*u.Mm,-1.27*u.Mm,1.27*u.Mm,0.0*u.Mm,25.5*u.Mm] #grid size


paper1 = {'photo_scale': 0.6*u.Mm,
          'chrom_scale': 0.1*u.Mm,
          'corona_scale': 2.5e3*u.Mm,         #scale height for the corona
          'coratio': 0.03*u.one,
          'model': 'paper1',
          'phratio': 0.0*u.one,
          'pixel': 0.36562475*u.Mm,              #(HMI pixel)
          'radial_scale': 0.10979002*u.Mm,
          'nftubes': 1,
          'B_corona': 9.2e-4*u.T,
          'pBplus': 1e-3*u.T}
paper1['chratio'] = 1*u.one - paper1['coratio'] - paper1['phratio']
paper1['Nxyz'] = [128,128,432] # 3D grid
paper1['xyz']  = [-1.27*u.Mm,1.27*u.Mm,-1.27*u.Mm,1.27*u.Mm,0.*u.Mm,8.62*u.Mm] #grid size
# uncomment to produce comaparable data to mfe_setup
#paper1['Nxyz'] = [127,128,128] # 3D grid
#paper1['xyz']  = [-1*u.Mm,1*u.Mm,-1*u.Mm,1*u.Mm,3.5e-2*u.Mm,1.6*u.Mm] #grid size

paper2a = {'photo_scale': 0.6*u.Mm,
           'chrom_scale': 0.1*u.Mm,
           'corona_scale': 2.5e3*u.Mm,         #scale height for the corona
           'coratio': 0.03*u.one,
           'model': 'paper2a',
           'phratio': 0.0*u.one,
           'pixel': 0.36562475*u.Mm,              #(HMI pixel)
           'radial_scale': 0.10979002*u.Mm,
           'nftubes': 4,
           'B_corona': 8.4e-4*u.T,
           'pBplus': 1e-3*u.T}
paper2a['chratio'] = 1*u.one - paper2a['coratio'] - paper2a['phratio']
paper2a['Nxyz'] = [160,80,432] # 3D grid
paper2a['xyz']  = [-1.59*u.Mm,1.59*u.Mm,-0.79*u.Mm,0.79*u.Mm,0.*u.Mm,8.62*u.Mm] #grid size

paper2b = {'photo_scale': 0.6*u.Mm,
           'chrom_scale': 0.1*u.Mm,
           'corona_scale': 2.5e3*u.Mm,         #scale height for the corona
           'coratio': 0.03*u.one,
           'model': 'paper2b',
           'phratio': 0.0*u.one,
           'pixel': 0.36562475*u.Mm,              #(HMI pixel)
           'radial_scale': 0.10979002*u.Mm,
           'nftubes': 4,
           'B_corona': 8.2e-4*u.T,
           'pBplus': 1.0e-3*u.T}
paper2b['chratio'] = 1*u.one - paper2b['coratio'] - paper2b['phratio']
paper2b['Nxyz'] = [50,50,140] # 3D grid
paper2b['xyz']  = [-0.49*u.Mm,0.49*u.Mm,-0.49*u.Mm,0.49*u.Mm,0*u.Mm,2.78*u.Mm] #grid size

paper2c = {'photo_scale': 0.6*u.Mm,
           'chrom_scale': 0.1*u.Mm,
           'corona_scale': 5e3*u.Mm,         #scale height for the corona
           'coratio': 0.03*u.one,
           'model': 'paper2c',
           'phratio': 0.0*u.one,
           'pixel': 0.36562475*u.Mm,              #(HMI pixel)
           'radial_scale': 0.10979002*u.Mm,
           'nftubes': 15,
           'B_corona': 5.95e-4*u.T,
           'pBplus': 1.0e-3*u.T}
paper2c['chratio'] = 1*u.one - paper2c['coratio'] - paper2c['phratio']
paper2c['Nxyz'] = [224,224,140] # 3D grid
paper2c['xyz']  = [-2.23*u.Mm,2.23*u.Mm,-2.23*u.Mm,2.23*u.Mm,0*u.Mm,2.78*u.Mm] #grid size

paper2d = {'photo_scale': 0.6*u.Mm,
           'chrom_scale': 0.1*u.Mm,
           'corona_scale': 5e3*u.Mm,         #scale height for the corona
           'coratio': 0.03*u.one,
           'model': 'paper2d',
           'phratio': 0.0*u.one,
           'pixel': 0.36562475*u.Mm,              #(HMI pixel)
           'radial_scale': 0.10979002*u.Mm,
           'nftubes': 18,
           'B_corona': 5.95e-4*u.T,
           'pBplus': 1.0e-3*u.T}
paper2d['chratio'] = 1*u.one - paper2d['coratio'] - paper2d['phratio']
paper2d['Nxyz'] = [224,224,140] # 3D grid
paper2d['xyz']  = [-2.23*u.Mm,2.23*u.Mm,-0.79*u.Mm,0.79*u.Mm,0*u.Mm,2.78*u.Mm] #grid size


def get_coords(Nxyz, xyz):
    """
    get_coords returns a non-dimensional dictionary describing the domain
    coordinates.
    """
    dz=(xyz[5]-xyz[4])/(Nxyz[2]-1)
    Z    = u.Quantity(np.linspace(xyz[4].value, xyz[5].value, Nxyz[2]), unit=xyz.unit)
    Zext = u.Quantity(np.linspace(Z.min().value-4.*dz.value, Z.max().value+4.*dz.value, Nxyz[2]+8), unit=Z.unit)
    coords = {
              'dx':(xyz[1]-xyz[0])/(Nxyz[0]-1),
              'dy':(xyz[3]-xyz[2])/(Nxyz[1]-1),
              'dz':(xyz[5]-xyz[4])/(Nxyz[2]-1),
              'xmin':xyz[0],
              'xmax':xyz[1],
              'ymin':xyz[2],
              'ymax':xyz[3],
              'zmin':xyz[4],
              'zmax':xyz[5],
              'Z':Z,
              'Zext':Zext
             }

    return coords
#-----------------------------------------------------------------------------

def get_hmi_coords(
                Nxyz,xyz,
                indx,
                dataset = 'hmi_m_45s_2014_07_06_00_00_45_tai_magnetogram_fits',
                sunpydir = os.path.expanduser(os.path.expanduser('~')+'/sunpy/data/'),
                figsdir = os.path.expanduser(os.path.expanduser('~')+'/figs/hmi/'),
                l_newdata = True,
                interpfactor=5,
                frac=[0,0],
                rank=0,
                lmpi=False
                   ):
    """
    get_coords returns a non-dimensional dictionary describing the domain
    coordinates.
    """
    dz=(xyz[5]-xyz[4])/(Nxyz[2]-1)
    Z    = u.Quantity(np.linspace(xyz[4].value,xyz[5].value,Nxyz[2]),
                      unit=u.Mm)
    Zext = u.Quantity(np.linspace(Z.min().value-4.*dz.value,
                                  Z.max().value+4.*dz.value, Nxyz[2]+8),
                      unit=u.Mm)
    s,x,y,FWHM,sdummy,xdummy,ydummy=get_hmi_map(
                indx,
                dataset = dataset,
                sunpydir = sunpydir,
                figsdir = figsdir,
                l_newdata = l_newdata,
                interpfactor=interpfactor,
                frac=frac,
                rank=rank,
                lmpi=lmpi
               )
    xmin=x.min()+(x.max()-x.min())*frac[0]
    xmax=x.min()+(x.max()-x.min())*(1-frac[0])
    ymin=y.min()+(y.max()-y.min())*frac[1]
    ymax=y.min()+(y.max()-y.min())*(1-frac[1])
    coords = {
              'dx':(xmax-xmin)/(Nxyz[0]-1),
              'dy':(ymax-ymin)/(Nxyz[1]-1),
              'dz':(xyz[5]-xyz[4])/(Nxyz[2]-1),
              'xmin':xmin.to(u.Mm),
              'xmax':xmax.to(u.Mm),
              'ymin':ymin.to(u.Mm),
              'ymax':ymax.to(u.Mm),
              'zmin':xyz[4],
              'zmax':xyz[5],
              'Z':Z,
              'Zext':Zext
             }

    return coords
#-----------------------------------------------------------------------------
#
def get_hmi_map(
        indx,
        dataset = 'hmi_m_45s_2014_07_06_00_00_45_tai_magnetogram.fits',
        sunpydir = os.path.expanduser(os.path.expanduser('~')+'/sunpy/data/'),
        figsdir = os.path.expanduser(os.path.expanduser('~')+'/figs/hmi/'),
        l_newdata = False,
        interpfactor=5,
        frac=[0,0],
        rank=0,
        lmpi=False
               ):
    """ indx is 4 integers: lower and upper indices each of x,y coordinates 
#    dataset of the form 'hmi_m_45s_2014_07_06_00_00_45_tai_magnetogram_fits'
#    """
    from scipy.interpolate import RectBivariateSpline
    from sunpy.net import vso
    import sunpy.map
    client = vso.VSOClient()
    results = client.query(vso.attrs.Time("2014/07/05 23:59:50",
                                          "2014/07/05 23:59:55"),
                           vso.attrs.Instrument('HMI'),
                           vso.attrs.Physobs('LOS_magnetic_field'))
    #print results.show()                       
    if lmpi:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD

    if l_newdata:
        if rank==0:
            if not os.path.exists(sunpydir):
                raise ValueError("in get_hmi_map set 'sunpy' dir for vso data\n"+
            "for large files you may want link to local drive rather than network"
            )
            client.get(results).wait(progress=True)
            lwait = False
        else:
            lwait = None
        if lmpi:
            lwait = comm.bcast(lwait, root=0)
 
    if not os.path.exists(figsdir):
        os.makedirs(figsdir)

    print sunpydir, dataset
    hmi_map = sunpy.map.Map(sunpydir+dataset)
    #hmi_map = hmi_map.rotate()
    #hmi_map.peek()
    s = hmi_map.data[indx[0]:indx[1],indx[2]:indx[3]] #units of Gauss Bz
    nx = s.shape[0]
    ny = s.shape[1]
    nx_int, ny_int = interpfactor*(nx-1)+1, interpfactor*(ny-1)+1 # size of interpolant 
    #pixel size in arc seconds
    if not l_astropy_0:
        dx, dy = hmi_map.scale.items()[0][1],hmi_map.scale.items()[1][1]
        x_int, y_int = np.mgrid[
                                hmi_map.xrange[0]+indx[0]*dx:
                                hmi_map.xrange[0]+(indx[1]-1)*dx:1j*nx_int,
                                hmi_map.yrange[0]+indx[2]*dy:
                                hmi_map.yrange[0]+(indx[3]-1)*dy:1j*ny_int
                               ]
        x, y = np.mgrid[
                        hmi_map.xrange[0]+indx[0]*dx:
                        hmi_map.xrange[0]+(indx[1]-1)*dx:1j*nx,
                        hmi_map.yrange[0]+indx[2]*dy:
                        hmi_map.yrange[0]+(indx[3]-1)*dy:1j*ny
                       ]
    else:
        dx, dy = hmi_map.scale.x.value,hmi_map.scale.y.value
        x_int, y_int = np.mgrid[
                                hmi_map.xrange[0].value+indx[0]*dx:
                                hmi_map.xrange[0].value+indx[1]*dx:1j*nx_int,
                                hmi_map.yrange[0].value+indx[2]*dy:
                                hmi_map.yrange[0].value+indx[3]*dy:1j*ny_int
                               ]
        x, y = np.mgrid[
                        hmi_map.xrange[0].value+indx[0]*dx:
                        hmi_map.xrange[0].value+indx[1]*dx:1j*nx,
                        hmi_map.yrange[0].value+indx[2]*dy:
                        hmi_map.yrange[0].value+indx[3]*dy:1j*ny
                       ]
        #arrays to interpolate s from/to
    fx = np.linspace(x_int.min(),x_int.max(),nx)
    fy = np.linspace(y_int.min(),y_int.max(),ny)
    xnew = np.linspace(x_int.min(),x_int.max(),nx_int)
    ynew = np.linspace(y_int.min(),y_int.max(),ny_int)
    f  = RectBivariateSpline(fx,fy,s)
    #The initial model assumes a relatively small region, so a linear
    #Cartesian map is applied here. Consideration may be required if larger
    #regions are of interest, where curvature or orientation near the lim
    #of the surface is significant. 
    s_int  = f(xnew,ynew) #interpolate s and convert units to Tesla
    interp_scale = 1. #contribution to field strength from location
    interp_scale += 2. #contribution to field strength from 4 lateral neighbours
    interp_scale += 0.146524 #contribution from nearest 4 diagonal neighbours
    xq = u.Quantity(x * 7.25e5, unit= u.m)
    yq = u.Quantity(y * 7.25e5, unit= u.m)
    xq_int = u.Quantity(x_int * 7.25e5, unit= u.m)
    yq_int = u.Quantity(y_int * 7.25e5, unit= u.m)
    sq = u.Quantity(s * 1e-4, unit= u.T)
    sq_int = u.Quantity(s_int * 1e-4 / interp_scale, unit= u.T)

    dx *= 7.25e5 * u.m
    dy *= 7.25e5 * u.m
    pixel = (dx+dy)*0.5
    FWHM  = 2*pixel/interpfactor

    import matplotlib.pyplot as plt
    xmin=xq.to(u.Mm).min().value+(xq.to(u.Mm).max().value-xq.to(u.Mm).min().value)*frac[0]
    xmax=xq.to(u.Mm).min().value+(xq.to(u.Mm).max().value-xq.to(u.Mm).min().value)*(1-frac[0])
    ymin=yq.to(u.Mm).min().value+(yq.to(u.Mm).max().value-yq.to(u.Mm).min().value)*frac[0]
    ymax=yq.to(u.Mm).min().value+(yq.to(u.Mm).max().value-yq.to(u.Mm).min().value)*(1-frac[0])

    cmax = max(-sq_int.min().value, sq_int.max().value)
    cmin=-cmax
    
    plt.figure()
    plt.pcolormesh(xq_int.T.to(u.Mm).value, yq_int.T.to(u.Mm).value, sq_int.T.value,
                   vmin=cmin, vmax=cmax)
    plt.xlabel('lon [Mm]')
    plt.ylabel('lat [Mm]')
    plt.axis([xmin,xmax,ymin,ymax])
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$B_z$ [T]')
    cbar.solids.set_edgecolor("face")
    plt.savefig('model_domain_hmi.png')
    plt.close()

    plt.figure()
    plt.pcolormesh(xq_int.T.to(u.Mm).value, yq_int.T.to(u.Mm).value, sq_int.T.value,
                   vmin=cmin, vmax=cmax)
    plt.axis([xq.to(u.Mm).value.min(),xq.to(u.Mm).value.max(),yq.to(u.Mm).value.min(),yq.to(u.Mm).value.max()])
    plt.xlabel('lon [Mm]')
    plt.ylabel('lat [Mm]')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$B_z$ [T]')
    cbar.solids.set_edgecolor("face")
    plt.savefig('interpolated_hmi.png')
    plt.close()

    plt.figure()
    plt.pcolormesh(xq.T.to(u.Mm).value, yq.T.to(u.Mm).value, sq.T.value,
                   vmin=cmin, vmax=cmax)
    plt.axis([xq.to(u.Mm).value.min(),xq.to(u.Mm).value.max(),yq.to(u.Mm).value.min(),yq.to(u.Mm).value.max()])
    plt.xlabel('lon [Mm]')
    plt.ylabel('lat [Mm]')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$B_z$ [T]')
    cbar.solids.set_edgecolor("face")
    plt.savefig('hmi.png')
    plt.close()

    return sq_int, xq_int, yq_int, FWHM, sq, xq, yq
