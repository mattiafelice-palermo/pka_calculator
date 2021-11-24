#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 17:09:55 2021

@author: mpalermo
"""
import os
from jobprocessor import JobProcessor
from pka_calculator import pKaCalculator

pka_jobs = [
    pKaCalculator(xyz.strip(".xyz"), 0, 0).calculate_pka
    for xyz in os.listdir("./xyz_files")
]

jobprocessor = JobProcessor(pka_jobs, maxcores=8, cores_per_job=2)
results = jobprocessor.run()
