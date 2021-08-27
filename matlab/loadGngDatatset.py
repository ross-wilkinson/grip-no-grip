#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 14:30:42 2021

@author: rosswilkinson

"""

import pandas as pd
from pymer4.models import Lmer

df = pd.read_csv('https://github.com/ross-wilkinson/datasets/blob/main/gripNoGrip.csv')
df = df.rename(columns={0: "Subject", 1: "Condition", 2: "Trial", 3: "MaxPowerCycle"})

model = Lmer('MaxPowerCycle ~ Condition + (1+Condition|Subject)', data=df)
mdf = model.fit()
model.summary()

