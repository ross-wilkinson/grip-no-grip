#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 22 02:11:33 2021

@author: rosswilkinson
"""
from pyomeca import Markers, Analogs
import numpy as np
import os

def c3d2trc(path,filename):
    
    # Get marker data from c3d file
    markers = Markers.from_c3d(path + '/' + filename + '.c3d', prefix_delimiter=":")
    analogs = Analogs.from_c3d(path + '/' + filename + '.c3d')
    
    # Set labels and trajectories to vars
    trajLabels = markers.channel.data
    trajData = markers.data
    
    # Transform data from UQ Qualisys XYZ to OpenSim XYZ
    # UQ X = OS X, UQ Y = OS -Z, UQ Z = OS Y
    trajData[[1, 2],:,:] = trajData[[2, 1],:,:]
    trajData[2,:,:] = -trajData[2,:,:]
    
    os.chdir(path)
    # Create new filename
    with open(filename + '.trc', 'w') as f:
        # Create header information
        f.write('PathFileType\t4\t(X/Y/Z)\t%s.trc\n' % filename)
        f.write('DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\t'+
                'OrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n')
        f.write('%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\n' 
                % (analogs.rate, markers.rate, markers.data.shape[2],
                   len(markers.channel), markers.units, markers.rate, 
                   markers.first_frame, markers.last_frame+1))
        # Create data frame
        header4 = 'Frame#\tTime\t'
        header5 = '\t\t'
        formatText = '%d\t%.4f\t'
        startFrame = markers.first_frame
        endFrame = markers.last_frame
        frameList = np.arange(startFrame,endFrame+1)
        time = markers.time
        dataOut = np.column_stack((frameList,time))
        
        # Use for loop to stack header and marker data
        n = 0
        for mm in range(len(markers.channel)):
            header4 = header4 + trajLabels[mm] + '\t\t\t'
            header5 = header5 + 'X%s\tY%s\tZ%s\t' % (mm,mm,mm) 
            formatText = formatText + '%f\t%f\t%f\t'
            dataOut = np.column_stack((dataOut,trajData[:3,mm,:].T))
            
        f.write(header4 + '\n')
        f.write(header5 + '\n')
        
        for i in range(len(dataOut)):
            f.write(formatText % tuple(dataOut[i,:]))
            f.write('\n')
        f.close()