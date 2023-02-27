#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 09:12:59 2022

@author: chrisyoung
"""

# Plot 3C waveforms for one station for specified event

import obspy
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth
from obspy.geodetics import kilometers2degrees
from obspy.taup import TauPyModel
import matplotlib.pyplot as plt

c = Client('IRIS')

# Get event info for time windowing
event_time = obspy.UTCDateTime('2011-03-11T05:46:23.2')

cat = c.get_events(starttime = event_time - 10, endtime = event_time + 10, 
                   minmagnitude = 9)

# Get station info for calculating predicted arrival time
inv = c.get_stations(network = 'TA', station = '637A', level='response')

st = c.get_waveforms(network='TA', station = '637A', location = '',
                     channel='BH?', starttime = event_time - 60, 
                     endtime = event_time + 3600)

# Get station coordinates
coords = inv.get_coordinates('TA.637A..BHE')

# Get source location
origin = cat[0].preferred_origin()

# Calculate epicentral distance to get predicted arrivals
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
st2 = st.slice(first_arrival - 60, first_arrival + 300)
# st2.plot(fig = fig)

# Deconvolve instrument response and plot windowed waveforms
st3 = st2.copy()
st3.attach_response(inv)      #attach instrument response info to wfm object
st3.remove_response(pre_filt = (1.0 / 10.0, 1.0 / 5.0, 1.0, 2.0), output = 'VEL')
fig = plt.figure()
st3.plot(fig = fig)

plt.figure(fig)
plt.savefig('Plot3C_Deconvolved.png', dpi = 300)


