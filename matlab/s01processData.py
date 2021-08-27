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

### Subject characteristics
number = 1
subName = 's01'
sex = 'Male'
age = 21 #yrs
mass = 81.5 #kg
height = 1.815 #m

### Get subject data file names
os.chdir(datPath)
subFileList = glob(subName + '*.c3d')
subFileList.remove('s01c0001.c3d')

# =============================================================================
# SCALE MODEL
# =============================================================================
### C3D 2 TRC
filename = 's01c0001'
osim.c3d2trc(datPath,filename)
### SETUP SCALE TOOL
osim.setupScaleTool(pathList, filename, subName, mass, height, age)
### -> RUN SCALE TOOL IN OPENSIM GUI

for files in subFileList:
    filename = files[:-4]

# =============================================================================
# SETUP INVERSE KINEMATICS
# =============================================================================
    ### C3D 2 TRC
    osim.c3d2trc(datPath,filename)
    ### SETUP IK TOOL
    osim.setupIkTool(datPath, modPath, resPath, setPath, name, filename)

# =============================================================================
# SETUP INVERSE DYNAMICS
# =============================================================================    
    ### WRITE EXTERNAL LOADS FILE
    osim.writeExternalLoadsFile(path, filename, rotateFlag, crankLength)  
    ### SETUP EXTERNAL LOADS
    osim.setupExternalLoads([pathList],filename)
    ### SETUP ID TOOL    
    osim.setupIdTool([pathList], filename)

# =============================================================================
# SETUP BODY KINEMATICS    
# =============================================================================
    ### SETUP ANALYZE TOOL: BODY KINEMATICS
    osim.setupBodyKinematics([pathList],filename)

# =============================================================================
# SETUP MOCO INVERSE
# =============================================================================  
    ### PROCESS EMG SIGNALS
    
    ### SETUP MOCO INVERSE
    



    