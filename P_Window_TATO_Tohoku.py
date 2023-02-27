#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 07:03:28 2023

@author: chrisyoung
"""

# Read TATO wfm data for Tohoku event, window P wave arrival and plot

import obspy
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
import matplotlib.pyplot as plt
from obspy.geodetics import gps2dist_azimuth
from obspy.geodetics import kilometers2degrees
from obspy.taup import TauPyModel

c = Client('IRIS')

# Set event time and fetch event info
event_time = obspy.UTCDateTime('2011-03-11T05:46:23.2')   #Tohoku event

cat = c.get_events(starttime = event_time - 10, endtime = event_time + 10, 
                    minmagnitude = 9)
inv = c.get_stations(network = 'IU', station = 'TATO', level='response')
starttime = UTCDateTime(2022,11,17,0,0,0)

# Fetch wfm data for 8 hours before event to 16 hours afterwards
st = c.get_waveforms(network='IU', station = 'TATO', location = '00',
                      channel='BHZ', starttime = event_time - (8 * 3600), 
                      endtime = event_time + 16 * 3600)

# Filter data and plot
tr = st[0];
tr.filter('bandpass', freqmin = 0.8, freqmax = 3.5)

# Get station coordinates
coords = inv.get_coordinates('IU.TATO..BHE')

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
print(first_arrival)

# Window around first arrival
tr2 = tr.slice(first_arrival - 60, first_arrival + 300)
fig = plt.figure()
tr2.plot(fig = fig)

plt.figure(fig)
plt.savefig('P_Window_TATO_Tohoku.png', dpi = 300)