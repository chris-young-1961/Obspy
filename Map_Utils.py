#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 17:37:52 2022

@author: chrisyoung
"""

# Package with map utilities translated from MatSeis

from math import sin, cos, asin, atan2, radians, degrees
from obspy.geodetics import gps2dist_azimuth
from obspy.geodetics import kilometers2degrees
import numpy as np

def great_circle(lat1, lon1, lat2, lon2, num_points):
    # Calculate dist, az between the 2 locations & convert to radians
    (dist,az,baz) = gps2dist_azimuth(lat1, lon1, lat2, lon2)
    dist_deg = kilometers2degrees(dist/1000.)
    dist_array = np.linspace(0.,radians(1) * dist_deg, num=num_points)
    az = radians(1) * az
    lat1 = radians(lat1)
    lon1 = radians(lon1)
    lat2 = radians(lat2)
    lon2 = radians(lon2)
    
    # Calculate trig.
    s1 = sin(lat1);
    s2 = sin(az);
    c1 = cos(lat1);
    c2 = cos(az);
    
    # Calculate points.
    latlon_list = []
    for dist in dist_array:
        s3 = sin(dist)
        c3 = cos(dist)
        lat = degrees(asin(s1*c3+c1*c2*s3))
        lon = degrees(lon1 + atan2(s2*s3,c1*c3-s1*c2*s3))
        latlon_list.append((lat, lon))
    return latlon_list

def small_circle(center_latitude, center_longitude, dist_degrees, num_points):
    az_array = radians(1) * (np.linspace(0.,360., num=num_points))
    dist = radians(dist_degrees)
    lat1 = radians(center_latitude)
    lon1 = radians(center_longitude)
    
    # Calculate trig.
    s1 = sin(lat1);
    s3 = sin(dist);
    c1 = cos(lat1);
    c3 = cos(dist);
    
    # Calculate points.
    latlon_list = []
    for az in az_array:
        s2 = sin(az)
        c2 = cos(az)
        lat = degrees(asin(s1*c3+c1*c2*s3))
        lon = degrees(lon1 + atan2(s2*s3,c1*c3-s1*c2*s3))
        latlon_list.append((lat, lon))
    return latlon_list