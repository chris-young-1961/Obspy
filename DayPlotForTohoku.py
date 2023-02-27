#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 09:12:59 2022

@author: chrisyoung
"""
# 1 day single station helicorder plot for 2011 Tohoku EQ 

import obspy
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
import matplotlib.pyplot as plt

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
fig = plt.figure()
tr.plot(fig = fig, type='dayplot',vertical_scaling_range = 0.02*max(abs(tr.data)))

plt.figure(fig)
plt.savefig('DayPlotForTohoku.png', dpi = 300)


