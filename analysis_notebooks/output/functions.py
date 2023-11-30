import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point

def blue2red_cmap(n, nowhite = False):
    """ combine two existing color maps to create a diverging color map with white in the middle
    n = the number of contour intervals
    """

    if (int(n/2) == n/2):
        # even number of contours
        nwhite=1
        nneg=n/2
        npos=n/2
    else:
        nwhite=2
        nneg = (n-1)/2
        npos = (n-1)/2

    if (nowhite):
        nwhite=0

    colors1 = plt.cm.Blues_r(np.linspace(0,1, int(nneg)))
    colors2 = plt.cm.YlOrRd(np.linspace(0,1, int(npos)))
    colorsw = np.ones((nwhite,4))

    colors = np.vstack((colors1, colorsw, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    return mymap



def contourmap_bothcontinents_robinson_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr,
 x1, x2, y1, y2, labels=True, fontsize=15, maskocean=False):
    """ plot a contour map of 2D data dat with coordinates lon and lat
        Input:
              fig = the figure identifier
              dat = the data to be plotted
              lon = the longitude coordinate
              lat = the latitude coordinate
              ci = the contour interval
              cmin = the minimum of the contour range
              cmax = the maximum of the contour range
              titlestr = the title of the map
              x1 = position of the left edge
              x2 = position of the right edge
              y1 = position of the bottom edge
              y2 = position of the top edge
              labels = True/False (ticks and  labels are plotted if true) 
              cmap = color map (only set up for blue2red at the moment)
    """

    # set up contour levels and color map
    nlevs = (cmax-cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)

    mymap = blue2red_cmap(nlevs)

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.Robinson(central_longitude=0))
    ax.set_aspect('auto')
    ax.add_feature(cfeature.COASTLINE, zorder=100)

    ax.set_title(titlestr, fontsize=fontsize)

    dat, lon = add_cyclic_point(dat, coord=lon)
    ax.contourf(lon, lat, dat, levels=clevs, cmap = mymap, extend="both", transform=ccrs.PlateCarree())

    ax.set_global()

    if (maskocean):
        ax.add_feature(cart.feature.OCEAN, zorder=1, edgecolor='k')

    return ax


def plotcolorbar(fig, ci, cmin, cmax, titlestr, x1, x2, y1, y2,
   orient='horizontal', posneg='both', ticks=None, fsize=14, nowhite=False,
   contourlines=False, contourlinescale=1):
    """plot a color bar
       Input:
           fig = the figure identified
           ci = the contour interval for the color map
           cmin = the minimum extent of the contour range
           cmax = the maximum extent of the contour range
           titlestr = the label for the color bar
           x1 = the location of the left edge of the color bar
           x2 = the location of the right edge of the color bar
           y1 = the location of the bottom edge of the color bar
           y2 = the location of the top edge of the color bar
           cmap = the color map to be used (only set up for blue2red at the moment)
           orient = the orientation (horizontal or vertical)
           posneg = if "both", both positive and negative sides are plotted
                    if "pos", only the positive side is plotted
                    if "net", only the negative side is plotted
           ticks = user specified ticklabels
           fsize = user specified font size
           contourlines = used to overplot contour lines
           contourlinescale = scale factor for contour lines to be overplotted
    """

    # set up contour levels and color map
    nlevs = (cmax-cmin)/ci + 1
    clevs = ci * np.arange(cmin/ci, (cmax+ci)/ci, 1)

    mymap = blue2red_cmap(nlevs, nowhite)

    clevplot=clevs
    if (posneg == "pos"):
        clevplot = clevs[clevs >= 0]
    if (posneg == "neg"):
        clevplot = clevs[clevs <= 0]

    norm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1])

    if (ticks):
        clb = mpl.colorbar.ColorbarBase(ax, cmap=mymap,
           orientation=orient, norm=norm, values=clevplot, ticks=ticks)
    else:
        clb = mpl.colorbar.ColorbarBase(ax, cmap=mymap,
           orientation=orient, norm=norm, values=clevplot)

    clb.ax.tick_params(labelsize=fsize)
    clb.set_label(titlestr, fontsize=fsize+2)

    if (contourlines):
        #clevlines = (clevs-ci/2.)*contourlinescale
        clevlines = clevs*contourlinescale
        clevlines = clevlines[np.abs(clevlines) > ci/2.]
        if (orient=='horizontal'):
            ax.vlines(clevlines[clevlines > 0],-5,5, colors='black', linestyle='solid')
            ax.vlines(clevlines[clevlines < 0],-5,5, colors='black', linestyle='dashed')
        if (orient=='vertical'):
            ax.hlines(clevlines[clevlines > 0],-10,15, colors='black', linestyle='solid')
            ax.hlines(clevlines[clevlines < 0],-10,15, colors='black', linestyle='dashed')


    return ax



def plotmap(fig, dat, ci, cmin, cmax, titlestr):
   
    cbartitlestr =  dat.attrs['units']


    ax = contourmap_bothcontinents_robinson_pos(fig, dat, dat.lon, dat.lat, ci, cmin, cmax, titlestr,
         0.05,0.45,0.75,0.97)

    ax = plotcolorbar(fig, ci, cmin, cmax,cbartitlestr,0.05,0.45,0.72,0.73)

    return ax

