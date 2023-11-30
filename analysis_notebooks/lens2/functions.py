import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import numpy as np
import cartopy as cart
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.ticker as mticker
from matplotlib.colors import ListedColormap ## used to create custom colormaps
import matplotlib.colors as mcolors
import matplotlib as mpl
from math import nan

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

def precip_cmap(n, nowhite=False):
    """ combine two existing color maps to create a diverging color map with white in the middle.
    browns for negative, blues for positive
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


    colors1 = plt.cm.YlOrBr_r(np.linspace(0,1, int(nneg)))
    colors2 = plt.cm.GnBu(np.linspace(0,1, int(npos)))
    colorsw = np.ones((nwhite,4))

    colors = np.vstack((colors1, colorsw, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    return mymap




def cosweightlonlat(darray,lon1,lon2,lat1,lat2, fliplon=True):
    """Calculate the weighted average for an [:,lat,lon] array over the region
    lon1 to lon2, and lat1 to lat2
    """
    # flip latitudes if they are decreasing
    if (darray.lat[0] > darray.lat[darray.lat.size -1]):
        print("flipping latitudes")
        darray = darray.sortby('lat')

    # flip longitudes if they start at -180
    if (fliplon):
        if (darray.lon[0] < 0):
            print("flipping longitudes")
            darray.coords['lon'] = (darray.coords['lon'] + 360) % 360
            darray = darray.sortby(darray.lon)


    region=darray.sel(lon=slice(lon1,lon2),lat=slice(lat1,lat2))
    weights = np.cos(np.deg2rad(region.lat))
    regionw = region.weighted(weights)
    regionm = regionw.mean(("lon","lat"))

    return regionm


def contourmap_bothcontinents_robinson_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr,
 x1, x2, y1, y2, labels=True, cmap="blue2red", fontsize=15, maskocean=False):
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

    if (cmap == "blue2red"):
        mymap = blue2red_cmap(nlevs)

    if (cmap == "precip"):
        mymap = precip_cmap(nlevs)

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


def contourmap_northamerica_fill_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr,
 x1, x2, y1, y2, labels=True, cmap="blue2red", maskocean=False, contourlines=False, contourlinescale=1):
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

    if (cmap == "blue2red"):
        mymap = blue2red_cmap(nlevs)

    if (cmap == "precip"):
        mymap = precip_cmap(nlevs)

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.PlateCarree())
    ax.set_aspect('auto')
    ax.add_feature(cfeature.COASTLINE)
    ax.set_extent([-170, -50, 10, 80], crs = ccrs.PlateCarree())

    if (labels):
        ax.set_xticks([-150, -100, -50], crs = ccrs.PlateCarree())
        ax.set_xticklabels(['150W','100W','50W'], fontsize=12)
        ax.set_yticks([20,40,60,80], crs = ccrs.PlateCarree())
        ax.set_yticklabels(['20N','40N','60N','80N'], fontsize=12)
        ax.xformatter = LongitudeFormatter()
        ax.yformatter = LatitudeFormatter()

    ax.set_title(titlestr, fontsize=16)

    dat, lon = add_cyclic_point(dat, coord=lon)
    ax.contourf(lon, lat, dat, levels=clevs, cmap = mymap, extend='both')

    if (maskocean):
        ax.add_feature(cart.feature.OCEAN, zorder=1, edgecolor='k')

    if (contourlines):
        clevlines = clevs*contourlinescale
        clevlines = clevlines[np.abs(clevlines) > ci/2. ]
        ax.contour(lon,lat,dat, levels=clevlines, colors='black', transform=ccrs.PlateCarree())


    return ax





def plotcolorbar(fig, ci, cmin, cmax, titlestr, x1, x2, y1, y2,
   cmap='blue2red', orient='horizontal', posneg='both', ticks=None, fsize=14, nowhite=False,
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

    if (cmap == "blue2red"):
        mymap = blue2red_cmap(nlevs, nowhite)

    if (cmap == "precip"):
        mymap = precip_cmap(nlevs, nowhite)

    if (cmap == "precip_nowhite"):
        mymap = precip_cmap_nowhite(nlevs)


    if (cmap == 'red2blue'):
        mymap = red2blue_cmap(nlevs, nowhite)

    clevplot=clevs
    if (posneg == "pos"):
        clevplot = clevs[clevs >= 0]
    if (posneg == "neg"):
        clevplot = clevs[clevs <= 0]

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1])
    norm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)

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





def plot_tmap_ensemblemean(fig, obs, model, start_year_of_decade):
    """ Plot the ensemble means """
    if (start_year_of_decade < 2030):
        ci = 0.2 ; cmax=3
    elif ((start_year_of_decade >= 2030) & (start_year_of_decade < 2060)):
        ci = 0.5 ; cmax=5
    elif ((start_year_of_decade >= 2060) & (start_year_of_decade < 2100)):
        ci = 1 ; cmax=10


        
    obsdat = obs.sel(year=slice(start_year_of_decade, start_year_of_decade+9)).mean('year')
    modeldat = model.sel(year=slice(start_year_of_decade, start_year_of_decade+9)).mean(['year','M'])

    obsbase = obs.sel(year=slice(1980,1989)).mean('year')
    modelbase = model.sel(year=slice(1980,1989)).mean(['year','M'])
   
    if ( (start_year_of_decade + 9) > np.max(obs.year.values) ) :
        obsdat[:] = nan

    if ( (start_year_of_decade + 9) > np.max(model.year.values) ) :
        obsdat[:] = nan
        
        

    ax = contourmap_bothcontinents_robinson_pos(fig, obsdat - obsbase,obsdat.lon,obsdat.lat,
        ci,-1*cmax,cmax,'Observations, '+str(start_year_of_decade)+'-'+str(start_year_of_decade + 9)+
         ' $-$ 1980-1989',0.05,0.45,0.75,0.97)

    ax = contourmap_bothcontinents_robinson_pos(fig, modeldat - modelbase,
            modeldat.lon, modeldat.lat, ci, -1*cmax, cmax, 'Model, '+
            str(start_year_of_decade)+'-'+str(start_year_of_decade + 9)+
            ' $-$ 1980-1989', 0.5,0.9,0.75,0.97)

    ax = plotcolorbar(fig, ci, -1*cmax, cmax, 'Temperature anomaly (K)',
               0.3,0.65,0.72,0.73)    

    return ax

def plot_prmap_ensemblemean(fig, obs, model, start_year_of_decade):
    """ Plot the ensemble means """
    obsdat = obs.sel(time=slice(start_year_of_decade, start_year_of_decade+9)).mean('time')
    modeldat = model.sel(time=slice(start_year_of_decade, start_year_of_decade+9)).mean(['time','M'])

    obsbase = obs.sel(time=slice(1980,1989)).mean('time')
    modelbase = model.sel(time=slice(1980,1989)).mean(['time','M'])
  
    if (start_year_of_decade < 2030): 
        ci = 0.1 ; cmax=1
    elif ((start_year_of_decade >= 2030) & (start_year_of_decade < 2060)):
        ci = 0.1 ; cmax=1
    elif ((start_year_of_decade >= 2060) & (start_year_of_decade < 2100)):
        ci = 0.1 ; cmax=1
    

    if ((start_year_of_decade + 9) > np.max(obs.time.values)):
        obsdat[:] = nan
    
 
    ax = contourmap_northamerica_fill_pos(fig, obsdat - obsbase,obsdat.lon,obsdat.lat,
        ci,-1*cmax,cmax,'Obs, '+str(start_year_of_decade)+'-'+str(start_year_of_decade + 9)+
         ' $-$ 1980-1989',0.05,0.25,0.75,0.95, cmap='precip')

    ax = contourmap_northamerica_fill_pos(fig, modeldat - modelbase,
            modeldat.lon, modeldat.lat, ci, -1*cmax, cmax, 'Model, '+
            str(start_year_of_decade)+'-'+str(start_year_of_decade + 9)+
            ' $-$ 1980-1989', 0.3,0.5,0.75,0.95, cmap='precip')

    ax = plotcolorbar(fig, ci, -1*cmax, cmax, 'Precipitation anomaly (mm/day)',
               0.05,0.5,0.72,0.73, cmap='precip')    

    return ax






def plot_tmap_members(fig, model, start_year_of_decade, members):
    """ plot individual ensemble members """

    x1 = [0.05,0.37,0.69,0.05,0.37,0.69,0.05,0.37,0.69]
    x2 = [0.32,0.64,0.95,0.32,0.64,0.95,0.32,0.64,0.95]
    y1 = [0.8,0.8,0.8,0.59,0.59,0.59,0.38,0.38,0.38]
    y2 = [0.95,0.95,0.95,0.74,0.74,0.74,0.53,0.53,0.53]


    if (start_year_of_decade < 2030): 
        ci = 0.2 ; cmax=3
    elif ((start_year_of_decade >= 2030) & (start_year_of_decade < 2060)):
        ci = 0.5 ; cmax=5
    elif ((start_year_of_decade >= 2060) & (start_year_of_decade < 2100)):
        ci = 1 ; cmax=10

    modelbase = model.sel(year=slice(1980,1989)).mean('year')
    modeldat = model.sel(year=slice(start_year_of_decade, start_year_of_decade + 9)).mean('year')

    for imem in np.arange(0,len(members)):
        ax = contourmap_bothcontinents_robinson_pos(fig,modeldat.sel(M=members[imem]) - modelbase.sel(M=members[imem]),
                 modeldat.lon, modeldat.lat, ci, -1*cmax, cmax, 'Member '+str(members[imem])+', '+
                                str(start_year_of_decade)+'-'+str(start_year_of_decade + 9),
                                x1[imem],x2[imem],y1[imem],y2[imem])


    ax = plotcolorbar(fig, ci, -1.*cmax, cmax, 'Temperature anomaly (K)',
            0.3,0.7,0.35,0.36) 


    return ax 


def plot_prmap_members(fig, model, start_year_of_decade, members):
    """ plot individual ensemble members """

    x1 = [0.05,0.28,0.51,0.74,
          0.05,0.28,0.51,0.74,
          0.05,0.28,0.51,0.74]
    x2 = [0.25,0.48,0.71,0.94,
          0.25,0.48,0.71,0.94,
          0.25,0.48,0.71,0.94]
    y1 = [0.8,0.8,0.8,0.8,
          0.55,0.55,0.55,0.55,
          0.3, 0.3,0.3,0.3]
    y2 = [1,1,1,1,
          0.75,0.75,0.75,0.75,
          0.5,0.5,0.5,0.5]


    if (start_year_of_decade < 2030): 
        ci = 0.1 ; cmax=1
    elif ((start_year_of_decade >= 2030) & (start_year_of_decade < 2060)):
        ci = 0.1 ; cmax=1
    elif ((start_year_of_decade >= 2060) & (start_year_of_decade < 2100)):
        ci = 0.1 ; cmax=1

    modelbase = model.sel(time=slice(1980,1989)).mean('time')
    modeldat = model.sel(time=slice(start_year_of_decade, start_year_of_decade + 9)).mean('time')

    for imem in np.arange(0,len(members)):
        ax = contourmap_northamerica_fill_pos(
                 fig,modeldat.sel(M=members[imem]) - modelbase.sel(M=members[imem]),
                 modeldat.lon, modeldat.lat, ci, -1*cmax, cmax, 'Member '+str(members[imem])+', '+
                                str(start_year_of_decade)+'-'+str(start_year_of_decade + 9),
                                x1[imem],x2[imem],y1[imem],y2[imem], cmap='precip')


    ax = plotcolorbar(fig, ci, -1.*cmax, cmax, 'Precipitation anomaly (mm/day)',
            0.3,0.7,0.27,0.28, cmap='precip') 


    return ax 





    


