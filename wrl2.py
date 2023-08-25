import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.visualization.wcsaxes import Quadrangle
from astropy.visualization.wcsaxes.frame import EllipticalFrame
from astropy.table import QTable
from astropy.coordinates import (SkyCoord, Distance, Galactic, 
                                 EarthLocation, AltAz)
from astropy.io import fits
from astropy.wcs import WCS

# Define sources

icrs_wcs = WCS({'naxis': 2,
                'naxis1': 324,
                'naxis2': 162,
                'crpix1': 162.5,
                'crpix2': 81.5,
                'cdelt1': -1,
                'cdelt2': 1,
                'ctype1': 'RA---AIT',
                'ctype2': 'DEC--AIT'})

def plot_sky_coordinates(skycoord,frame="icrs",**kwargs):
    
    c = skycoord.transform_to(frame)
    if frame=="galactic":
        x=c.l
        y=c.b
    elif frame=="icrs":
        x=c.ra
        y=c.dec
    else:
        raise ValueError("Only supported coordinate frames are 'icrs' and 'galactic'")
    plt.plot(-x.wrap_at(180 * u.deg).radian, y.radian, **kwargs)

for title, wcs in [('WR', icrs_wcs)]:
    plt.figure(figsize=(6, 4))

    ax = plt.subplot(111, projection=wcs, frame_class=EllipticalFrame)
    ax.set_title(title)

    # Manually set the plot limits to show the entire sphere
    ax.set_xlim(-0.5, 324-0.5)
    ax.set_ylim(-0.5, 162-0.5)
    ax.set_aspect(1.)

    # Native grid
    lon, lat = ax.coords
    lon.set_ticks(spacing=30*u.deg)
    lon.grid()
    lat.set_ticks(spacing=15*u.deg)
    lat.grid()

    pn = QTable.read('wr.csv')

    open_pn = SkyCoord(
        ra=pn['_RAJ2000'],
        dec=pn['_DEJ2000'],
        unit='deg')
    len(open_pn)

    ax.plot_coord(open_pn, '.')

    # Galactic plane
    galactic_plane = ax.get_coords_overlay('galactic')
    galactic_plane[1].set_ticks(spacing=180*u.deg)
    galactic_plane[1].grid(color='purple', ls='dashed')
    galactic_plane[0].set_ticks_visible(False)

    # Try to plot +/-60 deg in equatorial latitude
    spring1 = Quadrangle((326.9, 16.1)*u.deg, 33.1*u.deg, 21*u.deg, facecolor='blue', alpha=0.2,
                   transform=ax.get_transform('icrs'))
    ax.add_patch(spring1)

    spring2 = Quadrangle((0, 16.1)*u.deg, 42.9*u.deg, 27*u.deg, facecolor='blue', alpha=0.2,
                   transform=ax.get_transform('icrs'))
    ax.add_patch(spring2)

    fall = Quadrangle((108.1, 27.7)*u.deg, 71.9*u.deg, 11.9*u.deg, facecolor='green', alpha=0.2,
                   transform=ax.get_transform('icrs'))
    ax.add_patch(fall)

    fall2 = Quadrangle((-180, 27.7)*u.deg, 102.4*u.deg, 11.9*u.deg, facecolor='green', alpha=0.2,
                   transform=ax.get_transform('icrs'))
    ax.add_patch(fall2)

    fallb = Quadrangle((120, 39.6)*u.deg, 60*u.deg, 30.0*u.deg, facecolor='green', alpha=0.2,
                   transform=ax.get_transform('icrs'))
    ax.add_patch(fallb)

    fallb2 = Quadrangle((-180, 39.6)*u.deg, 75*u.deg, 30.0*u.deg, facecolor='green', alpha=0.2,
                   transform=ax.get_transform('icrs'))
    ax.add_patch(fallb2)

    fallc2 = Quadrangle((-105, 39.6)*u.deg, 27.4*u.deg, 7.5*u.deg, facecolor='green', alpha=0.2,
                   transform=ax.get_transform('icrs'))
    ax.add_patch(fallc2)

    l2gps1 = Quadrangle((32.0, 0.0)*u.deg, 182*u.deg, 5*u.deg, facecolor='red', alpha=0.2,
                   transform=ax.get_transform('galactic'))
    ax.add_patch(l2gps1)

    l2gps2 = Quadrangle((32.0, 0.0)*u.deg, 182*u.deg, -5*u.deg, facecolor='red', alpha=0.2,
                   transform=ax.get_transform('galactic'))
    ax.add_patch(l2gps2)


plt.legend()
plt.show()
