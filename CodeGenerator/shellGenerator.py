#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 24 11:08:41 2018

@author: nneuman
"""

#import numpy as np
import vrrGenerator
import ascsGenerator
import hrrGenerator



LMAX_global = 2;



LaMax = 2
LbMax = 2
LcMax = 2
LdMax = 2;
orderMax = LaMax+LbMax+LcMax+LdMax+1

FileName = 'OS' + hrrGenerator.stringType(LaMax) + hrrGenerator.stringType(LbMax) + hrrGenerator.stringType(LcMax) + hrrGenerator.stringType(LdMax)

OSfile = open(FileName,'w+')

for La in range(1,LaMax+1):
    for Lc in range(1):
        for order in range(orderMax-La-Lc):
            
            tmpout = vrrGenerator.doVRR_1(La,order,FileName)
            #OSfile.write(vrrGenerator.doVRR_1(La,order))
            #OSfile.write(' ')
            #print(' ')

    for Lc in range(1,LcMax+1):
        for order in range(orderMax-La-Lc):
            
            tmpout = ascsGenerator.doascs(La,Lc,order,FileName)
            #OSfile.write(ascsGenerator.doascs(La,Lc,order))
            #print(' ')    

for Lb in range(1,LbMax+1):
    for Ld in range(1,LdMax+1):
        hrrGenerator.dohrr2(LaMax,Lb,LcMax,Ld,2,0,FileName)
        #OSfile.write(hrrGenerator.dohrr2(LaMax,Lb,LcMax,Ld,2,0))

OSfile.close()
