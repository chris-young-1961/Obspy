#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 10:13:40 2022

@author: chrisyoung
"""

# Create and plot specified travel time curves using obspy.taup

import numpy as np
from obspy.taup import TauPyModel
import matplotlib.pyplot as plt

m = TauPyModel(model = 'ak135')
phase_list = ['P','S','PcP','PP','PKP']
color_list = ['red','blue','green','orange','purple','brown','olive','pink','cyan','gray']
dist_range = np.arange(0., 180, 0.5)

fig = plt.figure()

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
    plt.plot(dist_valid,time_valid,color=color,linestyle='-')
    
ax = fig.axes[0]
ax.set_xlabel('distance [deg]')
ax.set_ylabel('time [sec]')

plt.figure(fig)
plt.savefig('PlotTtCurves.png', dpi = 300)