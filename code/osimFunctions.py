#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 22 02:11:33 2021

@author: rosswilkinson
"""

### Import libraries
from pyomeca import Markers, Analogs
import os
import numpy as np
import pandas as pd
from glob import glob
import xml.etree.ElementTree as ET
from scipy import signal
import warnings
import math


# =============================================================================
# ROTATE MARKER DATA        
# =============================================================================

def rotateMarkerData(path, filename):
    
    ### Get marker data from c3d file
    markers = Markers.from_c3d(path + '/' + filename + '.c3d', prefix_delimiter=":")
    
    ### Set labels and trajectories to vars
    trajLabels = markers.channel.data
    trajData = markers.data
    
    ### Create bike reference frame using markers on base of ergometer
    B2 = np.squeeze(trajData[:3,'LBIKE'==trajLabels,:])
    B3 = np.squeeze(trajData[:3,'RBIKE'==trajLabels,:])
    
    ### Get XYZ distance between markers
    dist = np.array(B2 - B3).mean(axis=1)
    
    ### Create new Y-axis in bike reference frame.
    newY = dist / np.linalg.norm(dist)
    
    ### Set Z-Axis (vertical) same as lab
    labZ = np.array([0, 0, 1])
    
    ### New X-Axis = cross product of newY and labZ
    newX = np.cross(newY, labZ)
    
    ### New Z-Axis = cross product of newY and newX
    newZ = np.cross(newX, newY)
    
    ### Create 3D rotation matrix: Lab to Bike
    rM3D = np.vstack((newX, newY, newZ))
    
    ### Rotate marker data: Lab to Bike
    trajDataRotated = np.empty_like(trajData)
    
    for i in range(len(trajLabels)):        
        trajDataRotated[:3,i,:] = np.dot(rM3D,trajData[:3,i,:])
    
    return trajLabels, trajDataRotated, markers


# =============================================================================
# C3D 2 TRC
# =============================================================================

def c3d2trc(path, filename, rotateFlag):
    
    ### Call ROTATEMARKERDATA function
    if rotateFlag:
        trajLabels, trajData, markers = rotateMarkerData(path, filename)
    else:
        ### Get marker data from c3d file
        markers = Markers.from_c3d(path + '/' + filename + '.c3d', prefix_delimiter=":")
        
        ### Set labels and trajectories to vars
        trajLabels = markers.channel.data
        trajData = markers.data
    
    ### Get analog data from c3d file
    analogs = Analogs.from_c3d(path + '/' + filename + '.c3d')
    
    ### Transform data from UQ Qualisys XYZ to OpenSim XYZ
    # UQ X = OS X, UQ Y = OS -Z, UQ Z = OS Y
    trajData[[1, 2],:,:] = trajData[[2, 1],:,:]
    trajData[2,:,:] = -trajData[2,:,:]
    
    
    ### Create new file
    os.chdir(path)
    with open(filename + '.trc', 'w') as f:
        # Create header information
        f.write('PathFileType\t4\t(X/Y/Z)\t%s.trc\n' % filename)
        f.write('DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\t'+
                'OrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n')
        f.write('%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\n' 
                % (analogs.rate, markers.rate, markers.data.shape[2],
                   len(markers.channel), markers.units, markers.rate, 
                   markers.first_frame, markers.last_frame+1))
        ### Create data frame
        header4 = 'Frame#\tTime\t'
        header5 = '\t\t'
        formatText = '%d\t%.4f\t'
        startFrame = markers.first_frame
        endFrame = markers.last_frame
        frameList = np.arange(startFrame,endFrame+1)
        time = markers.time
        dataOut = np.column_stack((frameList,time))
        
        ### Use for loop to stack header and marker data
        for mm in range(len(trajLabels)):
            header4 = header4 + trajLabels[mm] + '\t\t\t'
            header5 = header5 + 'X%s\tY%s\tZ%s\t' % (mm+1,mm+1,mm+1) 
            formatText = formatText + '%f\t%f\t%f\t'
            dataOut = np.column_stack((dataOut,trajData[:3,mm,:].T))
            
        f.write(header4 + '\n')
        f.write(header5 + '\n')
        
        for i in range(len(dataOut)):
            f.write(formatText % tuple(dataOut[i,:]))
            f.write('\n')
        f.close()


# =============================================================================
# SETUP SCALE TOOL
# =============================================================================

def setupScaleTool(pathList, filename, subName, mass, height, age):
    
    ### Set directory paths
    datPath = pathList[0]
    modPath = pathList[1]
    resPath = pathList[2]
    setPath = pathList[3]
    
    ### Call C3D2TRC function
    rotateFlag = False
    c3d2trc(datPath, filename, rotateFlag)
        
    ### Load marker data to get start and end times
    markers = Markers.from_c3d(datPath + '/' + filename + '.c3d', prefix_delimiter=":")
    startTime = float(markers.time[0].values)
    endTime = float(markers.time[-1].values)
    
    ### Parse generic setup file
    tree = ET.parse(modPath + '/' + 'setupScale.xml')
    root = tree.getroot()
    
    ### Edit ScaleTool
    ScaleTool = root.find('ScaleTool')
    ScaleTool.attrib = {'name': subName}
    ScaleTool.find('mass').text = str(mass)
    ScaleTool.find('height').text = str(height)
    ScaleTool.find('age').text = str(age)
    
    ### Edit ScaleTool: GenericModelMaker (GMM)
    GMM = ScaleTool.find('GenericModelMaker')
    GMM.attrib = {'name': subName}
    GMM.find('model_file').text = glob(modPath + '/*.osim')[0]
    GMM.find('marker_set_file').text = modPath + '/markerSet.xml'
    
    ### Edit ScaleTool: ModelScaler (MS)
    MS = ScaleTool.find('ModelScaler')
    MS.attrib = {'name': subName}
    MS.find('scaling_order').text = 'measurements'
    MS.find('MeasurementSet').attrib = {'file': modPath + '/measurementSet.xml'}
    MS.find('marker_file').text = datPath + '/' + filename + '.trc'
    MS.find('time_range').text = str([startTime, endTime])
    MS.find('preserve_mass_distribution').text = 'true'
    
    ### Edit ScaleTool: MarkerPlacer (MP)
    MP = ScaleTool.find('MarkerPlacer')
    MP.attrib = {'name': subName}
    MP.find('IKTaskSet').attrib = {'file': modPath + '/tasksScale.xml'}
    MP.find('marker_file').text = datPath + '/' + filename + '.trc'
    MP.find('time_range').text = str([startTime, endTime])
    MP.find('output_model_file').text = resPath + '/' + subName + 'modelScaled.osim'
    MP.find('output_motion_file').text = resPath + '/' + filename + '.mot'
    MP.find('output_marker_file').text = resPath + '/' + subName + 'markersScaled.osim'
    MP.find('max_marker_movement').text = str(-1)
    
    ### Write XML file
    tree.write(setPath + '/' + subName + 'setupScale.xml')
    
    
# =============================================================================
# SETUP IK TOOL    
# =============================================================================

def setupIkTool(pathList, filename, subName):
    
    ### Set directory paths
    datPath = pathList[0]
    modPath = pathList[1]
    resPath = pathList[2]
    setPath = pathList[3]
    
    ### Load marker data to get start and end times
    markers = Markers.from_c3d(datPath + '/' + filename + '.c3d', prefix_delimiter=":")
    startTime = float(markers.time[0].values)
    endTime = float(markers.time[-1].values)
    
    ### Parse generic setup file
    tree = ET.parse(modPath + '/' + 'setupInverseKinematics.xml')
    root = tree.getroot()
    
    ### Edit InverseKinematicsTool
    IKTool = root.find('InverseKinematicsTool')
    IKTool.attrib = {'name': filename}
    IKTool.find('results_directory').text = resPath
    IKTool.find('model_file').text = resPath + '/' + subName + 'modelScaled.osim'
    IKTool.find('constraint_weight').text = 'Inf'
    IKTool.find('accuracy').text = '1e-05'
    IKTool.find('IKTaskSet').attrib = {'file': modPath + '/tasksInverseKinematics.xml'}
    IKTool.find('marker_file').text = datPath + '/' + filename + '.trc'
    IKTool.find('time_range').text = str([startTime,endTime])
    IKTool.find('report_errors').text = 'true'
    IKTool.find('output_motion_file').text = resPath + '/' + filename + 'inverseKinematics.mot'
    IKTool.find('report_marker_locations').text = 'false'
    
    ### Write XML file
    tree.write(setPath + '/' + filename + 'setupInverseKinematics.xml')


# =============================================================================
# PROCESS CRANK SIGNALS
# =============================================================================

def processCrankSignals(path, filename, crankLength):
    
    ### Load analog data from c3d file
    analogs = Analogs.from_c3d(path + '/' + filename + '.c3d')
    
    ### Set labels and channels to vars
    analLabels = analogs.channel.values
    analData = analogs.data
    
    rTrq = analData['RTan'==analLabels]
    
    if 'RRad' in analLabels:
        rRad = analData['RRad'==analLabels,:]
    else:
        rRad = analData['RRadial'==analLabels,:]
        
    lTrq = analData['LTan'==analLabels,:]
    lRad = analData['LRad'==analLabels,:]
    
    ### Offset force signals at 2.5mV baseline
    offset = 2.5
    rTrq = rTrq - offset
    rRad = rRad - offset
    lTrq = lTrq - offset
    lRad = lRad - offset
    
    ### Set signal conversion factors
    rangeV = 5      #voltage
    rangeNm = 1000  #torque Axis1 = 1000, Axis2 = 700
    rangeN = 5000   #force Axis1 = 5000, Axis2 = 4000
    
    mv2nm = rangeNm / rangeV
    mv2n = rangeN / rangeV
    
    ### Convert signal voltages to SI Units (Newton meters, Newtons, Radians)
    rTrq = rTrq * mv2nm
    rRad = rRad * mv2n
    lTrq = lTrq * mv2nm
    lRad = lRad * mv2n
    
    ### Convert crank torque to tangential force
    rTan = rTrq / crankLength
    lTan = lTrq / crankLength
    
    ### Design lowpass filter
    order = 1
    fs = analogs.rate
    nyq = 0.5 * fs
    cutoff = 12 / nyq
    b, a = signal.butter(order, cutoff, 'low')
    
    ### Low-pass filter force signals
    rTanFilt = np.squeeze(signal.filtfilt(b, a, rTan))
    rRadFilt = np.squeeze(signal.filtfilt(b, a, rRad))
    lTanFilt = np.squeeze(signal.filtfilt(b, a, lTan))
    lRadFilt = np.squeeze(signal.filtfilt(b, a, lRad))
    
    return rTanFilt, rRadFilt, lTanFilt, lRadFilt
 
    
# =============================================================================
# CALCULATE CRANK ANGLE    
# =============================================================================

def calculateCrankAngle(path, filename, rotateFlag):
    
    ### Call ROTATEMARKERDATA function
    trajLabels, trajData, markers = rotateMarkerData(path, filename)
    
    ### Transform data from UQ Qualisys XYZ to OpenSim XYZ
    # UQ X = OS X, UQ Y = OS -Z, UQ Z = OS Y
    trajData[[1, 2],:,:] = trajData[[2, 1],:,:]
    trajData[2,:,:] = -trajData[2,:,:]
    
    nFrames = trajData.shape[2]
    
    ### Set marker vars
    rtoe = np.squeeze(trajData[:3,'rtoe'==trajLabels,:])
    rmt5 = np.squeeze(trajData[:3,'rmt5'==trajLabels,:])
    rcal = np.squeeze(trajData[:3,'rcal'==trajLabels,:])
    ltoe = np.squeeze(trajData[:3,'ltoe'==trajLabels,:])
    lmt5 = np.squeeze(trajData[:3,'lmt5'==trajLabels,:])
    lcal = np.squeeze(trajData[:3,'lcal'==trajLabels,:])   
    
    ### Check if any right-foot markers are missing
    if any(x.size == 0 for x in [rtoe, rmt5, rcal]):
        warnings.warn('Missing right foot markers')
        rSpindle = rmt5[:3,:]
    else:
        ### Get vector from toe to heel
        rFoot2Heel = rcal[:2,:] - rtoe[:2,:]
        ### Create perpendicular vector
        rPerp = rFoot2Heel[[1,0],:]
        rPerp[0,:] = -rPerp[0,:]        
        ### Create scaling factor
        rScale = 25 / np.linalg.norm(rPerp,axis=0)  
        ### Create scaled vector from foot to pedal (normal = 25 mm)
        rFoot2Pedal = rPerp * rScale       
        ### Calculate spindle XY coordinate
        rSpindle = rmt5[:2,:] + rFoot2Pedal             
        ### Calculate pedal Z coordinates
        zR = np.vstack((rtoe[2,:], rmt5[2,:], rcal[2,:])).mean(axis=0)
        rSpindle = np.vstack((rSpindle, zR))
    
    ### Check if any left-foot markers are missing
    if any(x.size == 0 for x in [ltoe, lmt5, lcal]):
        warnings.warn('Missing left foot markers')
        lSpindle = lmt5[:3,:]
    else:
        ### Get vector from toe to heel
        lFoot2Heel = lcal[:2,:] - ltoe[:2,:]    
        ### Create perpendicular vector
        lPerp = lFoot2Heel[[1,0],:]
        lPerp[0,:] = -lPerp[0,:]
        ### Create scaling factor
        lScale = 25 / np.linalg.norm(lPerp,axis=0)
        ### Create scaled vector from foot to pedal (normal = 25 mm)
        lFoot2Pedal = lPerp * lScale    
        ### Calculate spindle XY coordinate
        lSpindle = lmt5[:2,:] + lFoot2Pedal
        ### Calculate pedal Z coordinates
        zL = np.vstack((ltoe[2,:], lmt5[2,:], lcal[2,:])).mean(axis=0)
        lSpindle = np.vstack((lSpindle, zL))
    
    ### Calculate right-crank angle in global and crank XYZ
    angleGlobalRight = []
    angleTdcRight = []
    for i in range(nFrames):
        # CCW in Global XYZ
        # Set points
        a = np.copy(lSpindle[:2,i])
        a[0] = a[0] + 1
        b = lSpindle[:2,i]
        c = rSpindle[:2,i]
        
        # Calculate vectors
        ba = a-b
        bc = c-b
        
        # Calculate angle (0 to pi)
        cosAngle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
        angle = np.arccos(cosAngle)
        
        # Change decreasing values to >pi
        if bc[1] < 0:
            angleGlobalRight.append((2*math.pi) - angle)
        else:
            angleGlobalRight.append(angle)        
        
        # CW in Crank XYZ
        # Set points
        a = np.copy(lSpindle[:2,i])
        a[1] = a[1] + 1
        b = lSpindle[:2,i]
        c = rSpindle[:2,i]
        
        # Calculate vectors
        ba = a-b
        bc = c-b
        
        # Calculate angle (0 to pi)
        cosAngle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
        angle = np.arccos(cosAngle)
        
        # Change decreasing values to >pi
        if bc[0] < 0:
            angleTdcRight.append((2*math.pi) - angle)
        else:
            angleTdcRight.append(angle)
    
    ### Convert to numpy arrays
    angleGlobalRight = np.array(angleGlobalRight)
    angleTdcRight = np.array(angleTdcRight)
    
    ### Calculate left-crank angle as anti-phase signal of right   
    angleGlobalLeft = np.copy(angleGlobalRight)
    k1 = angleGlobalLeft > math.pi
    k2 = angleGlobalLeft < math.pi
    angleGlobalLeft[k1] = angleGlobalLeft[k1] - math.pi
    angleGlobalLeft[k2] = angleGlobalLeft[k2] + math.pi
    
    return angleTdcRight, angleGlobalRight, angleGlobalLeft, rSpindle.T, lSpindle.T, markers


# =============================================================================
# ROTATE CRANK FORCES
# =============================================================================
    
def rotateCrankForces(path, filename, rotateFlag, crankLength):
    
    ### Call PROCESSCRANKSIGNALS function
    rTan, rRad, lTan, lRad = processCrankSignals(path, filename, crankLength)
    
    ### Call CALCULATECRANKANGLE function
    _, thetaR, thetaL, pFR, pFL, markers = calculateCrankAngle(path, filename, rotateFlag)
    
    ### Pre-allocate output vars
    rFR = []
    rFL = []
    
    ### Create 2D Rotation Matrices: Crank to Bike
    for i in range(len(rTan)):
        # Right
        rMR = [ [np.cos(thetaR[i]), -np.sin(thetaR[i])],
               [np.sin(thetaR[i]), np.cos(thetaR[i])] ]
        
        # Left                                  
        rML = [ [np.cos(thetaL[i]), -np.sin(thetaL[i])], 
               [np.sin(thetaL[i]), np.cos(thetaL[i])] ]          
        
        # Set X-axis in crank reference frame as negative radial force
        fXR = -rRad[i]
        fXL = -lRad[i]
        
        # % Set Y-axis in crank reference frame as negative tangential force
        fYR = -rTan[i]
        fYL = -lTan[i]
        
        # Multiply forces by 2D rotation matrix
        fXYR = [fXR, fYR]
        fXYRBike = np.dot(rMR, fXYR)
        fXYL = [fXL, fYL]
        fXYLBike = np.dot(rML, fXYL)
        
        # Set Z-coordinate to zero. No forces in Z-axis.
        fXYZRBike = np.append(fXYRBike, 0)
        fXYZLBike = np.append(fXYLBike, 0)
        
        # Switch to reaction force for OpenSim
        rFR.append(-fXYZRBike)
        rFL.append(-fXYZLBike)
    
    # Convert to numpy arrays
    rFR = np.array(rFR)
    rFL = np.array(rFL)
    
    return rFR, rFL, pFR, pFL, markers


# =============================================================================
# WRITE EXTERNAL LOADS FILE
# =============================================================================

def writeExternalLoadsFile(path, filename, rotateFlag, crankLength):
    
    ### Call ROTATECRANKFORCES function
    rFR, rFL, pFR, pFL, markers = rotateCrankForces(path, filename, rotateFlag, crankLength)
    
    ### Create time vector
    time = markers.time
    
    ### Create data array
    dataOut = np.column_stack((time, rFR, pFR, rFL, pFL))
    
    ### Replace NaN with zero
    dataOut = np.nan_to_num(dataOut)
    
    ### Create column names
    colNames = ['time\t', 'forceRightX\t', 'forceRightY\t', 'forceRightZ\t', 
                'pointRightX\t', 'pointRightY\t', 'pointRightZ\t',
                'forceLeftX\t', 'forceLeftY\t', 'forceLeftZ\t', 
                'pointLeftX\t', 'pointLeftY\t', 'pointLeftZ\n']
    
    ### Create new file
    with open(path + '/' + filename + '.mot', 'w') as f:
        # Create header information
        f.write('External Loads File\n')
        f.write('version=1\n')
        f.write('nRows=%d\n' % len(time))
        f.write('nColumns=%d\n' % len(colNames))        
        f.write('Range=%s-%s seconds\n' % (time[0], time[-1]))
        f.write('endheader\n')
        f.write(str(colNames))
        
        # Write data array into file
        formatText = '%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f'
        for i in range(len(dataOut)):
            f.write(formatText % tuple(dataOut[i,:]))
            f.write('\n')         
        f.close()


# =============================================================================
# SETUP EXTERNAL LOADS
# =============================================================================

def setupExternalLoads(pathList, filename):
    
    ### Set directory paths
    datPath = pathList[0]
    modPath = pathList[1]
    resPath = pathList[2]
    setPath = pathList[3]
    
    ### Parse generic setup file
    tree = ET.parse(modPath + '/' + 'setupExternalLoads.xml')
    root = tree.getroot()
    
    ### Edit ExternalLoads
    EL = root.find('ExternalLoads')
    EL.attrib = {'name': filename}
    
    EF1 = EL.find('objects')[0]
    EF1.find('isDisabled').text = 'false'
    EF1.find('applied_to_body').text = 'calcn_r'
    EF1.find('force_expressed_in_body').text = 'ground'
    EF1.find('point_expressed_in_body').text = 'ground'
    EF1.find('force_identifier').text = 'forceRight'
    EF1.find('point_identifier').text = 'pointRight'
        
    EF2 = EL.find('objects')[1]
    EF2.find('isDisabled').text = 'false'
    EF2.find('applied_to_body').text = 'calcn_l'
    EF2.find('force_expressed_in_body').text = 'ground'
    EF2.find('point_expressed_in_body').text = 'ground'
    EF2.find('force_identifier').text = 'forceLeft'
    EF2.find('point_identifier').text = 'pointLeft'

    EL.find('datafile').text = setPath + '/' + filename + 'externalLoads.mot'
    EL.find('external_loads_model_kinematics_file').text = resPath + '/' + filename + 'inverseKinematics.mot'
    EL.find('lowpass_cutoff_frequency_for_load_kinematics').text = str(12)
    
    ### Write XML file
    tree.write(setPath + '/' + filename + 'setupExternalLoads.xml')  


# =============================================================================
# SETUP ID TOOL
# =============================================================================

def setupIdTool(pathList, filename, subName):
    
    ### Set directory paths
    datPath = pathList[0]
    modPath = pathList[1]
    resPath = pathList[2]
    setPath = pathList[3]
    
    ### Load marker data to get start and end times
    markers = Markers.from_c3d(datPath + '/' + filename + '.c3d', prefix_delimiter=":")
    startTime = float(markers.time[0].values)
    endTime = float(markers.time[-1].values)
    
    ### Parse generic setup file
    tree = ET.parse(modPath + '/' + 'setupInverseDynamics.xml')
    root = tree.getroot() 
    
    ### Edit InverseDynamicsTool
    IDTool = root.find('InverseDynamicsTool')
    IDTool.attrib = {'name': filename}
    IDTool.find('results_directory').text = resPath
    IDTool.find('model_file').text = resPath + '/' + subName + 'modelScaled.osim'
    IDTool.find('time_range').text = str([startTime,endTime])
    IDTool.find('forces_to_exclude').text = 'actuators muscles'
    IDTool.find('external_loads_file').text = setPath + '/' + filename + 'setupExternalLoads.xml'
    IDTool.find('coordinates_file').text = resPath + '/' + filename + 'inverseKinematics.mot'
    IDTool.find('lowpass_cutoff_frequency_for_coordinates').text = str(12)
    IDTool.find('output_gen_force_file').text = filename + 'inverseDynamics.sto'
    IDTool.find('joints_to_report_body_forces').text = 'All'
    
    ### Write XML file
    tree.write(setPath + '/' + filename + 'setupInverseDynamics.xml')  
    

# =============================================================================
# SETUP ANALYZE TOOL: BODY KINEMATICS
# =============================================================================

def setupBodyKinematics(pathList, filename, subName):
    
    ### Set directory paths
    datPath = pathList[0]
    modPath = pathList[1]
    resPath = pathList[2]
    setPath = pathList[3]
    
    ### Load marker data to get start and end times
    markers = Markers.from_c3d(datPath + '/' + filename + '.c3d', prefix_delimiter=":")
    startTime = float(markers.time[0].values)
    endTime = float(markers.time[-1].values)
    
    ### Parse generic setup file
    tree = ET.parse(modPath + '/' + 'setupAnalyze.xml')
    root = tree.getroot()
    
    ### Edit AnalyzeTool
    AnalTool = root.find('AnalyzeTool')
    AnalTool.attrib = {'name': filename}
    AnalTool.find('model_file').text = resPath + '/' + subName + 'modelScaled.osim'
    AnalTool.find('replace_force_set').text = 'false'
    AnalTool.find('results_directory').text = resPath
    AnalTool.find('output_precision').text = str(8)
    AnalTool.find('initial_time').text = str(startTime)
    AnalTool.find('final_time').text = str(endTime)
    AnalTool.find('solve_for_equilibrium_for_auxiliary_states').text = 'false'
    AnalTool.find('maximum_number_of_integrator_steps').text = str(20000)
    AnalTool.find('maximum_integrator_step_size').text = str(1)
    AnalTool.find('minimum_integrator_step_size').text = str(1e-8)
    AnalTool.find('integrator_error_tolerance').text = str(1e-5)
    AnalTool.find('coordinates_file').text = resPath + '/' + filename + 'inverseKinematics.mot'
    AnalTool.find('lowpass_cutoff_frequency_for_coordinates').text = str(12)
    
    ### Edit BodyKinematics
    BK = AnalTool.find('AnalysisSet').find('objects')[2]
    BK.attrib = {'name': 'BodyKinematics'}
    BK.find('on').text = 'true'
    BK.find('start_time').text = str(startTime)
    BK.find('end_time').text = str(endTime)
    BK.find('step_interval').text = str(1)
    BK.find('in_degrees').text = 'true'
    BK.find('bodies').text = 'center_of_mass'
    BK.find('express_results_in_body_local_frame').text = 'false'
    
    ### Write XML file
    tree.write(setPath + '/' + filename + 'setupBodyKinematics.xml')


# =============================================================================
# SETUP ANALYZE TOOL: MUSCLE ANALYSIS
# =============================================================================
 


# =============================================================================
# SETUP ANALYZE TOOL: POINT KINEMATICS
# =============================================================================



# =============================================================================
# PROCESS EMG SIGNALS
# =============================================================================

def processEmgSignals(path, filename, muscleList):
    
    ### Load EMG data from c3d file
    emg = Analogs.from_c3d(path + '/' + filename + '.c3d', 
                           suffix_delimiter=".", usecols=muscleList)
    
    ### Design bandpass filter
    order = 1 # 2*n
    fs = emg.rate
    nyq = 0.5 * fs
    lowcut = 10 / nyq
    highcut = 99 / nyq
    b, a = signal.butter(order, (lowcut, highcut), 'band')

    ### Remove DC offset
    emgCentered = emg.data - emg.data.mean(axis=1, keepdims=True)

    ### Bandpass filter
    emgFilt = signal.filtfilt(b, a, emgCentered)    
    
    ### Convert to absolute values
    emgAbs = np.abs(emgFilt)
    
    ### Lowpass filter (smooth) absolute values
    order = 2 #2*n
    lowcut = 5 / nyq
    b, a = signal.butter(order, lowcut, 'low')
    
    emgFiltFilt = signal.filtfilt(b, a, emgAbs)
    
    cols = np.append('Time', muscleList)
    dataOut = np.vstack((emg.time, emgFiltFilt))
    
    emgProcessed = pd.DataFrame(dataOut.T, columns = cols)
    
    return emgProcessed
    
# =============================================================================
# SETUP MOCO INVERSE
# =============================================================================



# =============================================================================
# END    
# =============================================================================
        
        
  
        
  
    
  
    
        
        
        
        
        
        
        
        
        
        
