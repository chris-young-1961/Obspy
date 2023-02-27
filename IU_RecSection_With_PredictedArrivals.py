#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 07:38:30 2022

@author: chrisyoung
"""

# Plots record section for event (Tohoku) for specified network (IU) and shows
# predicted arrivals for user specified phases.

import obspy
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth
from obspy.geodetics import kilometers2degrees
from obspy.taup import TauPyModel
import matplotlib.pyplot as plt

# Create IRIS client to fetch data
c = Client('IRIS')

# Get event information
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

##filter data
freqmin = '0.8'
freqmax = '3.5'
st_filter = st.filter('bandpass', freqmin=float(freqmin), freqmax=float(freqmax))

# Calculate epicentral distance for each station
for tr in st_filter:
    sta_str = (tr.stats.network + '.' + tr.stats.station + 
               '.' + tr.stats.location + '.' + tr.stats.channel )  
    coords = inv.get_coordinates(sta_str)
    distance = gps2dist_azimuth(coords['latitude'], coords['longitude'],
                                         origin.latitude, origin.longitude)[0]
    distance_deg = kilometers2degrees(distance/1000)
    tr.stats.distance = 1000. * distance_deg
    
# Plot record section
fig = plt.figure()
st_filter.plot(type='section', recordlength=2400, orientation = 'horizontal',
        dist_degree=False, scale=2.0,
        time_down=True, linewidth=.25, grid_linewidth=.25, show=False, fig=fig)

# Fix y limits to be in degrees & add title
ax = fig.axes[0]
ax.set_ylabel('Offset [deg]')
ax.invert_yaxis()
ax.set_title(UTC_str + ', IU Network, ' + freqmin + '-' + freqmax + 
             ' Hz Bandpass')

# Plot predicted arrivals
markercolor = (205/255, 133/255, 63/255)
textcolor = (0/255, 0/255, 128/255)
m = TauPyModel(model = 'ak135')
phase_list = ['P','S','PcP','PKP']
num_phase = len(phase_list)
for tr in st_filter:
    distance_deg = tr.stats.distance/1000.
    ax.text(5.0, distance_deg, tr.stats.station, rotation=0.0,
            color=textcolor, va="bottom", ha="left", fontsize=4, zorder=10)
    for phase in phase_list:
        # sta_pha = ('sta=' + st_filter[ista].stats.station + 
        #            ', dist=' + str(distance_deg) + ', phase=' + phase_list[iphase])
        # print(sta_pha)
        arrival = m.get_travel_times(distance_in_degree=distance_deg, phase_list=[phase],
                                    source_depth_in_km = origin.depth/1000.)
        if len(arrival) > 0:
            ax.plot(arrival[0].time, distance_deg,marker='|',
                    markerfacecolor=markercolor, markeredgecolor = markercolor, 
                    markersize = 3, zorder=11)
plt.figure(fig)
plt.savefig('IU_RecSection_With_PredictedArrivals.png', dpi = 300)



