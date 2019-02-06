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
import helperFunctions as hF


LMAX_global = 2;



LaMax = 2
LbMax = 2
LcMax = 1
LdMax = 1;
orderMax = LaMax+LbMax+LcMax+LdMax+1

FileName = 'OS' + hrrGenerator.stringType(LaMax) + hrrGenerator.stringType(LbMax) + hrrGenerator.stringType(LcMax) + hrrGenerator.stringType(LdMax) + '.m'

OSfile = open(FileName,'w+')

OSfile.write('function gout' + ' = ' + FileName[0:-2] + '(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues)\n')

for order in range(orderMax):
    strSintegrals = 's_s_s_s_' + str(order) +' = gSSSSNValues(:,' + str(order+1) + ');\n'
    OSfile.write(strSintegrals)
    
#s_s_s_s_0 = gSSSSNValues(:,1);
#s_s_s_s_1 = gSSSSNValues(:,2);
#s_s_s_s_2 = gSSSSNValues(:,3);
#s_s_s_s_3 = gSSSSNValues(:,4);
#s_s_s_s_4 = gSSSSNValues(:,5);
#
OSfile.write('PAx = RPAValues(:,1);\n')
OSfile.write('PAy = RPAValues(:,2);\n')
OSfile.write('PAz = RPAValues(:,3);\n')

OSfile.write('QCx = RQCValues(:,1);\n')
OSfile.write('QCx = RQCValues(:,2);\n')
OSfile.write('QCx = RQCValues(:,3);\n')

OSfile.write('WPx = RWPValues(:,1);\n')
OSfile.write('WPx = RWPValues(:,2);\n')
OSfile.write('WPx = RWPValues(:,3);\n')

OSfile.write('WQx = RWQValues(:,1);\n')
OSfile.write('WQx = RWQValues(:,2);\n')
OSfile.write('WQx = RWQValues(:,3);\n')

OSfile.write('ABx = RABValues(:,1);\n')
OSfile.write('ABy = RABValues(:,2);\n')
OSfile.write('ABz = RABValues(:,3);\n')

OSfile.write('CDx = RCDValues(:,1);\n')
OSfile.write('CDy = RCDValues(:,2);\n')
OSfile.write('CDz = RCDValues(:,3);\n')

OSfile.write('oo2p = 0.5./pValues;\n')
OSfile.write('oo2q = 0.5./qValues;\n')
OSfile.write('qoppq = qValues./ppqValues;\n')
OSfile.write('oo2pq = 0.5./ppqValues;\n')
#PAy = RPAValues(:,2);
#PAz = RPAValues(:,3);
#
#QCx = RQCValues(:,1);
#QCy = RQCValues(:,2);
#QCz = RQCValues(:,3);
#
#WPx = RWPValues(:,1);
#WPy = RWPValues(:,2);
#WPz = RWPValues(:,3);
#
#WQx = RWQValues(:,1);
#WQy = RWQValues(:,2);
#WQz = RWQValues(:,3);
#
#oo2p = 0.5./pValues;
#oo2q = 0.5./qValues;
#qoppq = qValues./ppqValues;
#oo2pq = 0.5./ppqValues;

OSfile.close()


for Lc in range(1):
    for La in range(1,LaMax+LbMax+1):
    
        for order in range(orderMax-La-Lc):
            
            tmpout = vrrGenerator.doVRR_1(La,order,FileName)
            #OSfile.write(vrrGenerator.doVRR_1(La,order))
            #OSfile.write(' ')
            #print(' ')

for Lc in range(1,LcMax+LdMax+1):
    for La in range(2,LaMax+LbMax+1):
        for order in range(orderMax-La-Lc):
            
            tmpout = ascsGenerator.doascs(La,Lc,order,FileName)
        hrrGenerator.contract(La,0,Lc,0,0,FileName) #indexNbr is actually not needed    
            #OSfile.write(ascsGenerator.doascs(La,Lc,order))
            #print(' ')    

#for La in range(LaMax,LaMax+LbMax+1):
for Lb in range(0,LbMax):            
    for Ld in range(1,LdMax+1):
        
        #hrrGenerator.dohrr2(La,0,LcMax+LdMax-Ld,Ld,4,500,FileName)
        hrrGenerator.dohrr2(LaMax+LbMax-Lb,Lb,LcMax+LdMax-Ld,Ld,4,500,FileName)
        
for Lb in range(1,LbMax+1):
    hrrGenerator.dohrr2(LaMax+LbMax-Lb,Lb,LcMax,LdMax,2,500,FileName)
        #OSfile.write(hrrGenerator.dohrr2(LaMax,Lb,LcMax,Ld,2,0))

OSfile = open(FileName,'a+')
OSfile.write('end')
OSfile.close()
