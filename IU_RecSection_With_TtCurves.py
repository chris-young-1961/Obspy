#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 07:38:30 2022

@author: chrisyoung
"""

# Plot record section with travel time curves for specified event and station network.

import obspy
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth
from obspy.geodetics import kilometers2degrees
from obspy.taup import TauPyModel
import numpy as np
import matplotlib.pyplot as plt

# Create IRIS client to fetch data
c = Client('IRIS')

# Set event parameters and get event information
UTC_str = '2011-03-11T05:46:23.2'     #Tohoku event
event_time = obspy.UTCDateTime(UTC_str)
cat = c.get_events(starttime = event_time - 10, endtime = event_time + 10, 
                    minmagnitude = 9)
origin = cat[0].preferred_origin()

# Get IU network station info
inv = c.get_stations(network = 'IU', station = '*', location = '00',
                      channel='BHZ',level='channel')


# Event time to 1 hour afterwards, all loc 00 BHZ channels
st = c.get_waveforms(network='IU', station = '*', location = '00',
                      channel='BHZ', starttime = event_time, 
                      endtime = event_time + 1. * 3600)
num_sta = len(st)

# Filter data
freqmin = '0.8'
freqmax = '3.5'
st_filter = st.filter('bandpass', freqmin=float(freqmin), freqmax=float(freqmax))

##Calculate epicentral distance for each station
for tr in st_filter:
    sta_str = (tr.stats.network + '.' + tr.stats.station + 
               '.' + tr.stats.location + '.' + tr.stats.channel )  
    coords = inv.get_coordinates(sta_str)
    distance = gps2dist_azimuth(coords['latitude'], coords['longitude'],
                                         origin.latitude, origin.longitude)[0]
    distance_deg = kilometers2degrees(distance/1000)
    tr.stats.distance = 1000. * distance_deg
    
# Plot record section, i.e. disance sorted waveforms
fig = plt.figure()

st_filter.plot(type='section', recordlength=2400, orientation = 'horizontal',
        dist_degree=False, scale=2.0,
        time_down=True, linewidth=.25, grid_linewidth=.25, show=False, fig=fig)

# Label waveforms with STA
textcolor = (0/255, 0/255, 128/255)
ax = fig.axes[0]
for tr in st_filter:
    distance_deg = tr.stats.distance/1000.
    ax.text(5.0, distance_deg, tr.stats.station, rotation=0.0,
            color=textcolor, va="bottom", ha="left", fontsize=4, zorder=10)

# Fix y limits to be in degrees & add title
ax.set_ylabel('Offset [deg]')
ax.invert_yaxis()
ax.set_title(UTC_str + ', IU Network, ' + freqmin + '-' + freqmax + 
             ' Hz Bandpass')

# Add travel time curves
# First create Taupy object
m = TauPyModel(model = 'ak135')
  
# Specify phases and colors  for curves         
phase_list = ['P','S','PcP','PKP']
color_list = ['red','blue','green','orange','purple','brown','olive','pink','cyan','gray']
dist_range = np.arange(0., 180, 0.5)

# Now calculate and plot travel time curves
for postion, phase in enumerate(phase_list):
    dist_valid = np.array([0])
    time_valid = np.array([0])
    for dist in dist_range:
        arrival = m.get_travel_times(distance_in_degree=dist, phase_list=[phase],
                                      source_depth_in_km = 0.0)
        if len(arrival):
            dist_valid = np.concatenate((dist_valid,[dist]))
            time_valid = np.concatenate((time_valid,[arrival[0].time]))           
    dist_valid = dist_valid[1:]
    time_valid = time_valid[1:]
    color = color_list[postion]
    plt.plot(time_valid,dist_valid,color=color,linestyle='-',linewidth=0.5, 
             zorder=0)
    
plt.figure(fig)
# Save plot to file
plt.savefig('IU_RecSection_With_TtCurves.png', dpi = 300)



