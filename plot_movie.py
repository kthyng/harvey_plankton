'''
Make animation of data together.
'''

import pandas as pd
import tabs
import xarray as xr
import scs
import os
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cmocean.cm as cmo

land_10m = cartopy.feature.NaturalEarthFeature('physical', 'land', '10m',
                                        edgecolor='face',
                                        facecolor=cartopy.feature.COLORS['land'])
pc = cartopy.crs.PlateCarree()
merc = cartopy.crs.Mercator(central_longitude=-85.0)

tz = 'UTC'
dstart = '2017-9-20'
dend = '2017-10-2'

figname = 'tail'

os.makedirs('figures/%s' % figname, exist_ok=True)

# controlling min/max by hand because FGB data skews it and isn't in view
smin = 22  # min((df1[key].min(), df3[key].min()))
smax = 34  # max((df1[key].max(), df3[key].max()))
# currents arrows
cwidth = 0.004  # 0.002
cscale = 1500  # 2000
# wind arrows
wscale = 200
wwidth = 0.007
# hfradar
dd = 2
# scatter
s = 150


## Read in data ##

# Buoys
fname = 'data/buoys.csv'
if not os.path.exists(fname):
    buoynames = ['D', 'W', 'B', 'F', 'V', '42019', '42035']
    buoys = pd.DataFrame()
    for buoy in buoynames:
        # read in timezone and then localize to None
        dft = tabs.read(buoy, dstart, dend, resample=('30T', 0, 'instant'), tz=tz).tz_localize(None)
        buoys = buoys.join(dft, how='outer')
    buoys.index.rename('Dates [UTC]', inplace=True)
    buoys.to_csv(fname)
else:
    buoys = pd.read_csv(fname, parse_dates=True, index_col=0)

# buoy locations
loc = 'https://raw.githubusercontent.com/kthyng/tabswebsite/master/includes/buoys.csv'
blocs = pd.read_csv(loc, index_col=0)

# HF Radar
fname = 'data/hfradar.nc'
if not os.path.exists(fname):
    loc = 'http://hfrnet-tds.ucsd.edu/thredds/dodsC/HFR/USEGC/6km/hourly/RTV/HFRADAR_US_East_and_Gulf_Coast_6km_Resolution_Hourly_RTV_best.ncd'
    hf = xr.open_dataset(loc)
    hf[['u','v']].sel(time=slice(dstart,dend)).to_netcdf(fname)
else:
    hf = xr.open_dataset(fname)


# SCS
resample = '5T'
fname = 'data/scs.csv'
if not os.path.exists(fname):
    loc = '/Users/kthyng/Documents/data/HRR September 2017/Leg1/PS18_09_Leg1_DiMarco_SCS/Sea-Bird-Thermosalinograph-(converted-ASCII-data)_20170922-223421.Raw'
    ft = scs.read_file(loc).resample(resample, base=0, label='left').mean()

    loc = '/Users/kthyng/Documents/data/HRR September 2017/Leg2/SCS_Point Sur/PS18_09_Leg2_Whilden_SCS/Sea-Bird-Thermosalinograph-(converted-ASCII-data)_20170927-162751.Raw'
    dft = scs.read_file(loc).resample(resample, base=0, label='left').mean()
    ft = ft.append(dft)

    loc = '/Users/kthyng/Documents/data/HRR September 2017/Leg3/PS18_09_Leg3_Campbell_SCS/Sea-Bird-Thermosalinograph-(converted-ASCII-data)_20170929-172721.Raw'
    dft = scs.read_file(loc).resample(resample, base=0, label='left').mean()
    ft = ft.append(dft)
    ft.to_csv(fname)
else:
    ft = pd.read_csv(fname, parse_dates=True, index_col=0)

# ADCP
fname = 'data/adcp.nc'
if not os.path.exists(fname):
    loc = '/Users/kthyng/Documents/data_processed/adcp/ps1809l1/ps1809l1_postproc/wh300.unrotated.phasepos.amp/contour/wh300.nc'
    adcp1 = xr.open_dataset(loc)

    loc = '/Users/kthyng/Documents/data_processed/adcp/ps1809l3/ps1809l3_postproc/wh300.unrotated.phasepos.phaseshift/contour/wh300.nc'
    adcp3 = xr.open_dataset(loc)

    adcp = xr.concat([adcp1,adcp3],dim='time')

    # u1 = adcp1.u.isel(depth_cell=0).rolling(time=20, center=True).mean().data
    # v1 = adcp1.v.isel(depth_cell=0).rolling(time=20, center=True).mean().data

    adcp.to_netcdf(fname)
else:
    adcp = xr.open_dataset(fname)
adcp['u'].isel(depth_cell=0).data = adcp.u.isel(depth_cell=0).rolling(time=20, center=True).mean().data
adcp['v'].isel(depth_cell=0).data = adcp.v.isel(depth_cell=0).rolling(time=20, center=True).mean().data

# Sampling locations (and times)
loc = 'data/leg1_latlon.csv'
ll = pd.read_csv(loc, index_col=0)
ll['lat'] = ll['DegLat'] + ll['MinLat']/60 + ll['SecLat']/3600
ll['lon'] = -(ll['DegLon'] + ll['MinLon']/60 + ll['SecLon']/3600)


## Make Images ##

# time for plot
start = pd.Timestamp('2017-9-23 00:00')  # in UTC
# dt = pd.Timedelta('1 hour')
dt = pd.Timedelta('15 min')
# end = pd.Timestamp('2017-9-24 00:00')
end = pd.Timestamp('2017-10-1 16:00')

counter = 0
date = start

while date < end:

    date = start + counter*dt
    fname = 'figures/%s/%s.png' % (figname,date.strftime('%Y-%m-%dT%H%M'))
    if os.path.exists(fname):
        counter += 1
        continue

    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_axes([0.06, 0.01, 0.93, 0.95], projection=merc)
    ax.set_extent([-97.8, -93.4, 27, 29.8], pc)
    gl = ax.gridlines(linewidth=0.2, color='gray', alpha=0.5, linestyle='-', draw_labels=True)
    # the following two make the labels look like lat/lon format
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabels_bottom = False  # turn off labels where you don't want them
    gl.ylabels_right = False
    ax.add_feature(land_10m, facecolor='0.8')
    ax.set_facecolor('0.8')

    # Add date
    datestr = date.strftime('%b %d %H:%M, %Y')
    ax.text(0.19, 0.93, datestr, fontsize=16, transform=ax.transAxes)

    # Add data — find nearest regardless of samping frequency of data

    # Label TABS buoys
    lons = np.asarray([blocs.loc[name,'lon'] for name in buoynames])
    lats = np.asarray([blocs.loc[name,'lat'] for name in buoynames])
    for lon, lat, col in zip(lons, lats, buoynames):
        ax.text(lon + 0.02, lat, col, transform=pc, fontsize=10, color='0.5')

    # buoy salinity
    cols = buoys.columns[['Salinity' in col for col in buoys.columns]]
    colbuoys = [col.split(':')[0] for col in cols]
    lons = [blocs.loc[name,'lon'] for name in colbuoys]
    lats = [blocs.loc[name,'lat'] for name in colbuoys]
    values = buoys.reindex([date], method='nearest').loc[:,cols]

    mappable = ax.scatter(lons, lats, c=values.values[0], s=s, cmap=cmo.haline,
                          transform=pc, vmin=smin, vmax=smax)
    cax = fig.add_axes([0.2, 0.85, 0.4, 0.02])
    cb = fig.colorbar(mappable, cax=cax, orientation='horizontal', pad=0.02, shrink=0.5, aspect=40, extend='both')
    cb.set_label('Salinity', fontsize=14)

    # buoy currents
    cols = buoys.columns[['East [cm/s]' in col for col in buoys.columns]]
    cols2 = buoys.columns[['North [cm/s]' in col for col in buoys.columns]]
    colbuoys = [col.split(':')[0] for col in cols]
    lons = np.asarray([blocs.loc[name,'lon'] for name in colbuoys])
    lats = np.asarray([blocs.loc[name,'lat'] for name in colbuoys])
    values = buoys.reindex([date], method='nearest').loc[:,cols]
    values2 = buoys.reindex([date], method='nearest').loc[:,cols2]

    Q = ax.quiver(lons, lats, values.values[0], values2.values[0], transform=pc,
                  width=cwidth, scale=cscale, headlength=2, headaxislength=2, zorder=10)
    qk = ax.quiverkey(Q, 0.075, 0.75, 50, '0.5 m$\cdot$s$^{-1}$ current', labelpos='E', coordinates='axes')


    # buoy wind
    cols = buoys.columns[['East [m/s]' in col for col in buoys.columns]]
    cols2 = buoys.columns[['North [m/s]' in col for col in buoys.columns]]
    colbuoys = [col.split(':')[0] for col in cols]
    lons = np.asarray([blocs.loc[name,'lon'] for name in colbuoys])
    lats = np.asarray([blocs.loc[name,'lat'] for name in colbuoys])
    values = buoys.reindex([date], method='nearest').loc[:,cols]
    values2 = buoys.reindex([date], method='nearest').loc[:,cols2]

    Q = ax.quiver(lons, lats, values.values[0], values2.values[0], transform=pc,
                  width=wwidth, scale=wscale, headlength=2, headaxislength=2, zorder=10,
                  alpha=1, color='0.6', linewidths=(0.7,), edgecolors=('k'))
    qk = ax.quiverkey(Q, 0.075, 0.63, 5, '5 m$\cdot$s$^{-1}$ wind', labelpos='E', coordinates='axes')


    # HF radar
    lat = hf['lat'].data[::dd]
    lon = hf['lon'].data[::dd]
    u = hf['u'].sel(time=date, method='nearest').data[::dd,::dd]*100  # cm/s
    v = hf['v'].sel(time=date, method='nearest').data[::dd,::dd]*100  # cm/s

    # matches the currents arrows above
    Q = ax.quiver(lon, lat, u, v, transform=pc,
                  width=cwidth, scale=cscale, headlength=2, headaxislength=2,
                  zorder=10, alpha=0.25)
    qk = ax.quiverkey(Q, 0.075, 0.7, 50, '0.5 m$\cdot$s$^{-1}$ current,\nHF Radar', labelpos='E', coordinates='axes')


    # SCS salt

    if figname == 'tail':  # keep showing scs for a time period
        dates = pd.date_range(start=date-pd.Timedelta('1 day'), end=date, freq=dt)
        ftt = ft.reindex(dates, method='nearest')
    else:
        ftt = ft.reindex([date], method='nearest')


    ax.scatter(ftt['lon'], ftt['lat'], c=ftt['Practical salinity'], s=s,
               marker='s', cmap=cmo.haline, zorder=1,
               transform=pc, vmin=smin, vmax=smax)


    # ADCP
    # if figname == 'tail':
    #     lat = adcp['lat'].sel(time=slice(date-pd.Timedelta('3 hours'),date)).data
    #     lon = adcp['lon'].sel(time=slice(date-pd.Timedelta('3 hours'),date)).data
    #     u = adcp['u'].sel(time=slice(date-pd.Timedelta('3 hours'),date)).isel(depth_cell=0).data*100  # cm/s
    #     v = adcp['v'].sel(time=slice(date-pd.Timedelta('3 hours'),date)).isel(depth_cell=0).data*100  # cm/s
    # else:
    lat = adcp['lat'].sel(time=date, method='nearest').data.reshape(1)
    lon = adcp['lon'].sel(time=date, method='nearest').data.reshape(1)
    u = adcp['u'].sel(time=date, method='nearest').isel(depth_cell=0).data.reshape(1)*100  # cm/s
    v = adcp['v'].sel(time=date, method='nearest').isel(depth_cell=0).data.reshape(1)*100  # cm/s

    # matches the currents arrows above
    Q = ax.quiver(lon, lat, u, v, transform=pc,
                  width=cwidth, scale=cscale, headlength=2, headaxislength=2, zorder=10)


    # label sampling station locations
    for station in ll.index:
        lon = ll['lon'][station]
        lat = ll['lat'][station]
        ax.plot(lon, lat, '.', transform=pc, zorder=0, color='0.5', markersize=1)
        ax.text(lon + 0.01, lat, station, transform=pc, fontsize=8, color='0.5')

    # add legend
    ax.plot(0.3, 0.76, 'o', markersize=10, color='#3E929A', transform=ax.transAxes)
    ax.text(0.33, 0.75, 'Buoy', color='k', transform=ax.transAxes)
    ax.plot(0.3, 0.71, 's', markersize=10, color='#3E929A', transform=ax.transAxes)
    ax.text(0.33, 0.7, 'Flow through data', color='k', transform=ax.transAxes)
    ax.plot(0.3, 0.66, '.', markersize=2, color='0.5', transform=ax.transAxes)
    ax.text(0.33, 0.65, 'Sampling station', color='k', transform=ax.transAxes)


    fig.savefig(fname, bbox_inches='tight', dpi=120)
    plt.close(fig)

    counter += 1
