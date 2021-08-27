#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 22 18:11:45 2021

@author: rosswilkinson
"""

# =============================================================================
# WRITE SETUP-SCALE FILE
# =============================================================================

import xml.etree.ElementTree as ET
import os
from glob import glob

def editScaleTool(datDir, modDir, resDir, setDir, name, mass, height, age, 
                  filename, startTime, endTime):
    ### Parse generic setup file
    os.chdir(modDir)
    tree = ET.parse('setupScale.xml')
    root = tree.getroot()
    
    ### Edit ScaleTool
    ScaleTool = root.find('ScaleTool')
    ScaleTool.attrib = {'name': name}
    ScaleTool.find('mass').text = mass
    ScaleTool.find('height').text = height
    ScaleTool.find('age').text = age
    
    ### Edit ScaleTool: GenericModelMaker (GMM)
    GMM = ScaleTool.find('GenericModelMaker')
    GMM.attrib = {'name': name}
    GMM.find('model_file').text = glob(modDir + '/*.osim')[0]
    GMM.find('marker_set_file').text = modDir + '/markerSet.xml'
    
    ### Edit ScaleTool: ModelScaler (MS)
    MS = ScaleTool.find('ModelScaler')
    MS.attrib = {'name': name}
    MS.find('scaling_order') = 'measurements'
    MS.find('MeasurementSet').attrib = {'file': modDir + '/measurementSet.xml'}
    MS.find('marker_file').text = datDir + '/' + filename + '.trc'
    MS.find('time_range').text = str([startTime, endTime])
    MS.find('preserve_mass_distribution') = 'true'
    
    ### Edit ScaleTool: MarkerPlacer (MP)
    MP = ScaleTool.find('MarkerPlacer')
    MP.attrib = {'name': name}
    MP.find('marker_file') = ''
    
    