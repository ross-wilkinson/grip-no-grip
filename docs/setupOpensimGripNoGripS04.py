#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 21 21:41:52 2021

@author: rosswilkinson
"""

# =============================================================================
# SETUP OPENSIM PIPELINE: GRIP NO GRIP: S01
# =============================================================================

### Set project folder paths
expPath = '/Users/rosswilkinson/Google Drive/projects/grip-no-grip'
datPath = expPath + '/data'
codPath = expPath + '/matlab'
modPath = expPath + '/model'
resPath = expPath + '/results'
setPath = expPath + '/setup'
docPath = expPath + '/docs'

pathList = [datPath, modPath, resPath, setPath]

### Import libraries
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from pyomeca import Markers, Analogs
from glob import glob
os.chdir(codPath)
import osimFunctions as osim
from scipy import signal
import xml.etree.ElementTree as ET

### Load subject data
os.chdir(docPath)
subjectData = pd.read_excel('subjectData.xlsx')

### Subject characteristics
subject = 4
if subject < 10:
    subName = 's0' + str(subject)
else:
    subName = 's' + str(subject)       
    
idx = subject - 1
sex = subjectData['sex'][idx]
age = subjectData['age'][idx] #yrs
mass = subjectData['mass'][idx] #kg
height = subjectData['height'][idx] / 100 #m

### Get subject data file names
os.chdir(datPath)
subFileList = glob(subName + '*.c3d')
subFileList.remove(subName + 'c0001.c3d')

# =============================================================================
# SCALE MODEL
# =============================================================================
### C3D 2 TRC
filename = subName + 'c0001'
rotateFlag = 0
osim.c3d2trc(datPath,filename,rotateFlag)
### SETUP SCALE TOOL
osim.setupScaleTool(pathList, filename, subName, mass, height, age)
### -> RUN SCALE TOOL IN OPENSIM GUI

for files in subFileList:
    filename = files[:-4]
    
# =============================================================================
# SETUP INVERSE KINEMATICS
# =============================================================================
    rotateFlag = 1   
    ### C3D 2 TRC
    osim.c3d2trc(datPath,filename,rotateFlag)
    ### SETUP IK TOOL
    osim.setupIkTool(pathList, filename, subName)

# =============================================================================
# SETUP INVERSE DYNAMICS
# =============================================================================    
    crankLength = 0.175
    ### WRITE EXTERNAL LOADS FILE
    osim.writeExternalLoadsFile(datPath, filename, rotateFlag, crankLength)  
    ### SETUP EXTERNAL LOADS
    osim.setupExternalLoads(pathList,filename)
    ### SETUP ID TOOL    
    osim.setupIdTool(pathList, filename, subName)

# =============================================================================
# SETUP BODY KINEMATICS    
# =============================================================================
    ### SETUP ANALYZE TOOL: BODY KINEMATICS
    osim.setupBodyKinematics(pathList, filename, subName)

# =============================================================================
# PROCESS EMG SIGNALS
# =============================================================================  
    ### PROCESS EMG SIGNALS
    muscleList = ['Biceps', 'Latissi Dorsi', 'Glut max', 'Erector spinal', 'VL', 'RF', 'BF', 'Gastroc med']
    emgProcessed = osim.processEmgSignals(datPath, filename, muscleList)
    os.chdir(resPath)
    emgProcessed.to_csv(filename + 'emgProcessed.csv', index=False)

    