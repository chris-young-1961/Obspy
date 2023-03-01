#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 15:30:03 2022

@author: chrisyoung
"""

# Plot waveform and spectragram for specified event and station

import obspy
from obspy.clients.fdsn import Client
import matplotlib.pyplot as plt
import numpy as np
from obspy.geodetics import gps2dist_azimuth
from obspy.geodetics import kilometers2degrees
from obspy.taup import TauPyModel
import matplotlib as mpl

c = Client('IRIS')

# Waveform window parameters
network = 'IU'
station = 'MAJO'
channel = 'BHZ'
lead_time = 5 * 60.   #seconds
total_time = 15 * 60.  #seconds

# Set event time and fetch event info. In this case event if found by searching
# both by time and position (radius around a specified point). Can also constrain
# by magnitude
event_str = 'NK Nuclear Explosion'
UTC_str = '2017-09-03T03:30:01.940'     #2017 NK Test
event_lat = 41.343
event_lon = 129.036
radius_deg = 180.
event_time = obspy.UTCDateTime(UTC_str)
cat = c.get_events(starttime = event_time - 10, endtime = event_time + 10, 
                    latitude=event_lat, longitude=event_lon, maxradius=0.1)
origin = cat[0].preferred_origin()    #only expecting one event in catalog

# Get metadata for the specified station.
inv = c.get_stations(network = 'IU', station = station, level='response')

# Get station coordinates
# coords = inv.get_coordinates('IU.MAJO..BHZ')
coords = inv.get_coordinates(network + '.' + station + '..' + channel)

# Get source location
origin = cat[0].preferred_origin()

# Calculate distance to get predicted arrivals
distance = gps2dist_azimuth(origin.latitude, origin.longitude,
                            coords['latitude'], coords['longitude'])
distance = kilometers2degrees(distance[0]/1000)

# Calculate predicted arrivals
m = TauPyModel(model = 'ak135')

arrivals = m.get_ray_paths(distance_in_degree=distance, 
                            source_depth_in_km = origin.depth/1000.)

#arrivals.plot()    #plot ray paths on Earth cross-section

first_arrival = origin.time + arrivals[0].time   #arrivals sorted by time

# Fetch wfm data for specified window around first arrival
st = c.get_waveforms(network=network, station = station, location = '00',
                      channel=channel, starttime = first_arrival - lead_time, 
                      endtime = first_arrival + (total_time - lead_time))
tr = st[0]
t = np.arange(0, tr.stats.npts / tr.stats.sampling_rate, tr.stats.delta)

# Create figure object
fig = plt.figure()

# Plot spectragram in top subplot
cmap = mpl.cm.cool
ax_top = fig.add_subplot(2,1,1)
tr.spectrogram(axes=ax_top, dbscale = True, log=True,cmap='jet', clip = [.5,1.0])
ax_top.set_xlim(min(t),max(t))
ax_top.set_ylabel('frequency [Hz]')
ax_top.set_title((event_str + ', ' + UTC_str + ', sta = ' + station +
'.' + channel),fontsize = 9)

# Plot waveform in bottom subplot
ax_bot = fig.add_subplot(2,1,2)
ax_bot.plot(t,tr.data,'k-', linewidth=0.5,)
ax_bot.set_xlim(min(t),max(t))
ax_bot.set_xlabel('seconds')

plt.figure(fig)
# Save figure to a file
plt.savefig('WaveformAndSpectragram.png', dpi = 300)