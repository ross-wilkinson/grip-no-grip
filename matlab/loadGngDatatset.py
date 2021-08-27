#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 14:30:42 2021

@author: rosswilkinson
"""

import pandas as pd

path = '/Users/rosswilkinson/Google Drive/projects/grip-no-grip/results'
filename = 'df.csv'

df = pd.read_csv(path + '/' + filename, header=None)
df = df.rename(columns={0: "Subject", 1: "Condition", 2: "Trial", 3: "MaxPowerCycle"})

