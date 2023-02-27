#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 07:38:30 2022

@author: chrisyoung
"""

# Plot station map and record setion for Feb. 2023 Turkey event.
# Can specify network to plot and station distance range.

import obspy
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth
from obspy.geodetics import kilometers2degrees
from obspy.taup import TauPyModel
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
# Next is a custom package that needs to be in same directory/folder
from Map_Utils import small_circle

# Function to strip STA out of NET.STA.LOC.CHAN
def get_sta(stachan):
    sta = []
    index = 3   #skip 2 character NET code and '.'
    while stachan[index] != '.':
        index += 1
    sta = stachan[3:(index)]
    return sta

# Set paramters for record section
freqmin = '1.0'                 # BP filter min frequency
freqmax = '4.0'                 # BP filter max frequency
min_dist_deg = 00.
max_dist_deg = 10.  
y_lim = [min_dist_deg,max_dist_deg]         # max station distance
time_duration = 1.00 * 3600     # time duration in seconds to plot after event time
lead_time = 5 * 60.             #lead time in seconds before event
scale_fac = (y_lim[1] - y_lim[0]) /30       # wfm amplitude scaling
font_size = 5                              # sta label fontsize
# Set parameters for map
small_circle_radii_deg = [5,10]          #r adii of small circles around event


# Create IRIS client to fetch data
c_iris = Client('IRIS')
c_koeri = Client('KOERI')
# network = 'IU'
# channel = 'BHZ'
# location = '00'
network = 'HL'
network = '*'
channel = 'HHZ,BHZ'
location = '*'
# Next is in case there are stations with bad data you don't want to plot
# Figure out the bad stations by plotting and seeing which waveforms don't look good.
exclude_stations = ['GAZ','KMRS','DIKM','MLAZ','DAT','CTKS','CTYL','EZN',
                    'ARMT','GAZK','ERIK','KCTX','KLYT','EDRB','GELI','TAHT',
                    'URFA','YOZ','KVT','BCK','VRTB','MDNY','TASB','SILT']

## read event information
# UTC_str = '2011-03-11T05:46:23.2'     #Tohoku event
# event_time = obspy.UTCDateTime(UTC_str)
# cat = c.get_events(starttime = event_time - 10, endtime = event_time + 10, 
#                     minmagnitude = 9)
# UTC_str = '2017-09-03T03:30:01.940'     #2017 NK Test
# event_lat = 41.343
# event_lon = 129.036
# radius_deg = 180.
UTC_str = '2023-02-06T01:17:35.000'     #2023 Turkey EQ
event_lat = 37.174
event_lon = 37.032

event_time = obspy.UTCDateTime(UTC_str)
cat = c_iris.get_events(starttime = event_time - 10, endtime = event_time + 10, 
                    latitude=event_lat, longitude=event_lon, maxradius=0.1)
origin = cat[0].preferred_origin()    #only expecting one event in catalog

# Get network station info for event time.
# Constraining time is necessary or you may get multiple items per station
inv = c_koeri.get_stations(network = network, station = '*', location = location,
                      channel=channel,level='channel', 
                      latitude = event_lat, longitude = event_lon,
                      minradius=min_dist_deg, maxradius=max_dist_deg,
                      starttime = event_time - 10, endtime = event_time + 10)

# Calculate epicentral distances to event
stachan_list = inv.get_contents()['channels']   #returns NET.STA.LOC.CHAN
distance_list = [];
coords_list = []
for sta_chan in stachan_list:
    coords = inv.get_coordinates(sta_chan)
    distance = gps2dist_azimuth(coords['latitude'], coords['longitude'],
                                         origin.latitude, origin.longitude)[0]
    distance_deg = kilometers2degrees(distance/1000.)
    distance_list.append(distance_deg)
    coords_list.append((coords['longitude'], coords['latitude']))
sorted_list = sorted(zip(distance_list,stachan_list,coords_list))

# Find stations at distance <= max_dist_deg
counter = 0;
for dist_chan in sorted_list:
    if dist_chan[0] > max_dist_deg:
        break;
    else:
        counter += 1
        
trimmed_list = sorted_list[0:counter]
      
starttime = event_time - lead_time
endtime = starttime + time_duration
# Form string of stations (not NET.STA.LOC.CHAN)
sta_str = get_sta(trimmed_list[0][1])    #load first station, then skip this one in for loop
for stachan in trimmed_list[1:]:
    sta_str = sta_str + ',' + get_sta(stachan[1])
        
# Get waveforms from event time to 1 hour afterwards, all loc 00 BHZ channels.
# Unfortunately get_waveforms does not return traces in station order specified
st = c_koeri.get_waveforms(network=network, station = sta_str, location = location,
                      channel=channel, starttime = starttime, 
                      endtime = endtime)
st.resample(40)          #downsample

# Because get_waveforms doesn't keep station order, we need figure out 
# distances for each trace for record section plotting
for tr in st:
    # stachan_str = (tr.stats.network + '.' + tr.stats.station + 
    #                '.' + tr.stats.location + '.' + tr.stats.channel ) 
    # print(stachan_str)
    for dist_chan in trimmed_list:
        if tr.stats.station == get_sta(dist_chan[1]):
            tr.stats.distance = dist_chan[0]
            tr.stats.longitude = dist_chan[2][0]
            tr.stats.latitude = dist_chan[2][1]
            break

# Demean, taper and filter data
st_filter = st.detrend(type='demean')
st_filter = st_filter.taper(type='hann',max_percentage=0.05)
st_filter = st_filter.filter('bandpass', freqmin=float(freqmin), freqmax=float(freqmax))

# Setup record section plot
textcolor = (0/255, 0/255, 128/255)
fig = plt.figure(figsize=[10, 5])
ax1 = fig.add_subplot(1, 2, 2)   #right side

# Plot travel time curves 
m = TauPyModel(model = 'ak135')         
#phase_list = ['P','S','PcP','PP','PKP']
phase_list = ['P','S']
color_list = ['red','blue','green','orange','purple','brown','olive','pink','cyan','gray']
dist_range = np.arange(min_dist_deg, max_dist_deg, (max_dist_deg - min_dist_deg)/50.)
for postion, phase in enumerate(phase_list):
    dist_valid = np.array([0])
    time_valid = np.array([0])
    for dist in dist_range:
        arrival = m.get_travel_times(distance_in_degree=dist, phase_list=[phase],
                                      source_depth_in_km = 0.0)
        if len(arrival):      #if we got a TT prediction for this phase at this distance
            dist_valid = np.concatenate((dist_valid,[dist]))
            time_valid = np.concatenate((time_valid,[arrival[0].time]))           
    dist_valid = dist_valid[1:]       #skip the first item
    time_valid = time_valid[1:]
    color = color_list[postion]
    ax1.plot(lead_time + time_valid,dist_valid,color=color,linestyle='-',linewidth=0.5, 
             zorder=0)    #low zorder so everthing else gets drawn on top

## add waveforms
for tr in st_filter:
    if tr.stats.station in exclude_stations:
        continue
    else:
        data = tr.stats.distance + scale_fac * tr.data/(max(abs(tr.data)))
        time = (tr.stats.starttime - starttime + 
                np.linspace(0,tr.stats.endtime-tr.stats.starttime+tr.stats.delta,len(data)))
        ax1.plot(time,data,'k-', linewidth=0.2, figure=fig)
    
# Set axes limits and add labels
ax1.set_ylabel('distance [deg]')
ax1.set_xlabel('time [sec]')
ax1.set_xlim(0, time_duration)
ax1.set_ylim(min_dist_deg,max_dist_deg)
ax1.invert_yaxis()
ax1.set_title(UTC_str + ', ' + network + ' Net, ' + freqmin + '-' + freqmax + 
             ' Hz BP',fontsize=8)

# Label waveforms with STA
textcolor = (200/255, 0/255, 0/255)
for tr in st_filter:
    if tr.stats.station in exclude_stations:
        continue
    else:
        ax1.text(5.0, tr.stats.distance, 
                 tr.stats.station + '.' + tr.stats.channel, rotation=0.0,
                 color=textcolor, va="bottom", ha="left", fontsize=font_size, zorder=10)
    
# Add map, left side
min_latitude = 28.
max_latitude = 50.
min_longitude = 22.
max_longitude = 48.
ax2 = fig.add_subplot(1, 2, 1, projection=ccrs.PlateCarree())
ax2.set_extent([min_longitude, max_longitude, min_latitude, max_latitude], crs=ccrs.PlateCarree())
ax2.add_feature(cfeature.OCEAN, zorder=0)
ax2.add_feature(cfeature.LAND, zorder=0, edgecolor='black',linewidth=0.2)
# ax2.set_global()
ax2.gridlines(color=(0.5,0.5,0.5))     #lat lon lines

# Plot small circles around event (if any)
if any(small_circle_radii_deg):
    for sc_deg in small_circle_radii_deg:
        sc_pts = small_circle(origin.latitude, origin.longitude, sc_deg, 50)
        sc_lat, sc_lon = zip(*sc_pts)
        ax2.plot(sc_lon, sc_lat,color='k',linestyle='--',linewidth=0.5, 
                 transform=ccrs.PlateCarree(),zorder=5)

# Plot event
event_color = (220/255, 20/255, 60/255)   #crimson
ax2.plot(origin.longitude, origin.latitude,
          marker='*', markersize=6, markerfacecolor=event_color, 
          markeredgecolor = 'k', markeredgewidth=0.1, 
          transform=ccrs.PlateCarree(),
          zorder=10)

# Plot stations
sta_color = (255/255, 140/255, 0/255)   #dark orange
for tr in st_filter:
    if tr.stats.station in exclude_stations:
        continue
    else:
    # Plot station
        ax2.plot(tr.stats.longitude, tr.stats.latitude,
                  marker='^', markersize=4, markerfacecolor=sta_color, 
                  markeredgecolor = 'k',markeredgewidth=0.1,
                  transform=ccrs.PlateCarree(), zorder=5)
        # # Path from event to station
        # ax2.plot([origin.longitude, tr.stats.longitude],
        #           [origin.latitude, tr.stats.latitude], 
        #           color='black', linewidth=0.5,linestyle='--',
        #           transform=ccrs.PlateCarree(),zorder = 6)
        # Label station
        ax2.text(tr.stats.longitude, tr.stats.latitude, 
                  tr.stats.station,
                  va="bottom", ha="left", fontsize=5,
                  transform=ccrs.PlateCarree(),zorder = 7)  
    
plt.figure(fig)
plt.savefig('EventRecSectionPlusMap_Turkey2023_KOERI.png', dpi = 300)