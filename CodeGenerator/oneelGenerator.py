#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 24 09:55:15 2018

@author: nneuman

This routine generates code for one-electron (overlap and kinetic energy) integrals
Overlap recursion relation
[a+1i|b] = PAi*[a|b]+ai*(1/2p)*[a-1i|b]+bi*(1/2p)*[a|b-1i]
Kinetic Energy Recursion Relation
[a+1i|T|b] = PAi*[a|T|b]+ai*(1/2p)*[a-1i|T|b]+bi*(1/2p)*[a|T|b-1i]+2*(a*b/p)*([a+1i|b]-ai*(1/2p)*[a-1i|b])
"""

import numpy as np
import helperFunctions as hF


def doOneEl(La,Lb,order,indexNbr,FileName):
    OSfile = open(FileName,'a+')
    #I require something that reverts to vrr if Lc == 0
    shell_a = hF.shell(La)
    shell_b = hF.shell(Lb)
#    print(shell_a)
    
    for ta in range(np.size(shell_a,0)):
        current_a = shell_a[ta]
        for tb in range(np.size(shell_b,0)):
            current_b = shell_b[tb]
            
            if indexNbr == 1:
                increaseDir = hF.findIncreaseDir(current_a)
        
    #            print('current_a = ' + str(current_a))        
    #            print('increasedir = ' + str(increaseDir))
    #            print(hF.decrementDir(current_a,increaseDir,1))
    #            print('current_b = ' + str(current_b))       
        
        #This part contains the actual recursion relations (i.e. the logic)--------------------------
                lhsOverlap = np.array([current_a,current_b,[order]])
                PATermOverlap = np.array([hF.decrementDir(current_a,increaseDir,1),current_b,[order]])
                #WQTerm = np.array([current_a,hF.decrementDir(current_b,increaseDir,1),[order+1]])
        #        print(decrementDir(current,increaseDir,2))
                try:
                    am2TermOverlap = np.array([hF.decrementDir(current_a,increaseDir,2),current_b,[order]])
                except:
                    am2TermOverlap = np.array([[hF.decrementDir(current_a,increaseDir,2)],current_b,[order]])
    
                try:
                    am1bm1TermOverlap = np.array([hF.decrementDir(current_a,increaseDir,1),hF.decrementDir(current_b,increaseDir,1),[order]])
                    
                except:
                    try:
                        am1bm1TermOverlap = np.array([hF.decrementDir(current_a,increaseDir,1),[hF.decrementDir(current_b,increaseDir,1)],[order]])
                    except:
                        am1bm1TermOverlap = np.array([[hF.decrementDir(current_a,increaseDir,1)],hF.decrementDir(current_b,increaseDir,1),[order]])    
                
        #-------------------------------------------------------------------------------------------        
                u = ('x','y','z')
                
                lhTerm = 'S' + hF.ind2str(lhsOverlap) + ' = '
                firstTerm = 'PA' + u[increaseDir] + '*S'  + hF.ind2str(PATermOverlap)
                    
                if current_a[increaseDir]-1 > 0:
                    secondTermOverlap = ' + ' + str(current_a[increaseDir]-1) + '*oo2p*(S' + hF.ind2str(am2TermOverlap) + ')'
                else:
                    secondTermOverlap = ''
                if current_b[increaseDir] > 0: #Not the same as for index_a
                    thirdTermOverlap = ' + ' + str(current_b[increaseDir]) + '*oo2p*(S' + hF.ind2str(am1bm1TermOverlap) + ')'
                else:
                    thirdTermOverlap = ''
                stringOutput = lhTerm + firstTerm + secondTermOverlap + thirdTermOverlap + ';\n'
                OSfile.write(stringOutput)
                
            elif indexNbr == 2:
                increaseDir = hF.findIncreaseDir(current_b)
        
    #            print('current_a = ' + str(current_a))        
    #            print('increasedir = ' + str(increaseDir))
    #            print(hF.decrementDir(current_a,increaseDir,1))
    #            print('current_b = ' + str(current_b))       
        
        #This part contains the actual recursion relations (i.e. the logic)--------------------------
                lhsOverlap = np.array([current_a,current_b,[order]])
                PBTermOverlap = np.array([current_a,hF.decrementDir(current_b,increaseDir,1),[order]])
                #WQTerm = np.array([current_a,hF.decrementDir(current_b,increaseDir,1),[order+1]])
        #        print(decrementDir(current,increaseDir,2))
                try:
                    bm2TermOverlap = np.array([current_a,hF.decrementDir(current_b,increaseDir,2),[order]])
                except:
                    bm2TermOverlap = np.array([current_a,[hF.decrementDir(current_b,increaseDir,2)],[order]])
    
                try:
                    am1bm1TermOverlap = np.array([hF.decrementDir(current_a,increaseDir,1),hF.decrementDir(current_b,increaseDir,1),[order]])
                    
                except:
                    try:
                        am1bm1TermOverlap = np.array([hF.decrementDir(current_a,increaseDir,1),[hF.decrementDir(current_b,increaseDir,1)],[order]])
                    except:
                        am1bm1TermOverlap = np.array([[hF.decrementDir(current_a,increaseDir,1)],hF.decrementDir(current_b,increaseDir,1),[order]])    
                
        #-------------------------------------------------------------------------------------------        
                u = ('x','y','z')
                
                lhTerm = 'S' + hF.ind2str(lhsOverlap) + ' = '
                firstTerm = 'PB' + u[increaseDir] + '*S'  + hF.ind2str(PBTermOverlap)
                    
                if current_b[increaseDir]-1 > 0:
                    secondTermOverlap = ' + ' + str(current_b[increaseDir]-1) + '*oo2p*(S' + hF.ind2str(bm2TermOverlap) + ')'
                else:
                    secondTermOverlap = ''
                if current_a[increaseDir] > 0: #Not the same as for index_a
                    thirdTermOverlap = ' + ' + str(current_a[increaseDir]) + '*oo2p*(S' + hF.ind2str(am1bm1TermOverlap) + ')'
                else:
                    thirdTermOverlap = ''
                stringOutput = lhTerm + firstTerm + secondTermOverlap + thirdTermOverlap + ';\n'
                OSfile.write(stringOutput)                
            #print(stringOutput)
    #        print('lhs = '+ str(lhs))
    #        print('PATerm = '+ str(PATerm))
    #        print('WPTerm = '+ str(WPTerm))
    
    #Now we do the kinetic energy part
    for ta in range(np.size(shell_a,0)):
        current_a = shell_a[ta]
        for tb in range(np.size(shell_b,0)):
            current_b = shell_b[tb]
            
            if indexNbr == 1:
                increaseDir = hF.findIncreaseDir(current_a)
        
    #            print('current_a = ' + str(current_a))        
    #            print('increasedir = ' + str(increaseDir))
    #            print(hF.decrementDir(current_a,increaseDir,1))
    #            print('current_b = ' + str(current_b))       
        
        #This part contains the actual recursion relations (i.e. the logic)--------------------------
                lhsKinetic = np.array([current_a,current_b,[order]])
                PATermKinetic = np.array([hF.decrementDir(current_a,increaseDir,1),current_b,[order]])
                TermOverlap = np.array([current_a,current_b,[order]])
                #WQTerm = np.array([current_a,hF.decrementDir(current_b,increaseDir,1),[order+1]])
        #        print(decrementDir(current,increaseDir,2))
                try:
                    am2TermKinetic = np.array([hF.decrementDir(current_a,increaseDir,2),current_b,[order]])
                except:
                    am2TermKinetic = np.array([[hF.decrementDir(current_a,increaseDir,2)],current_b,[order]])
    
                try:
                    am1bm1TermKinetic = np.array([hF.decrementDir(current_a,increaseDir,1),hF.decrementDir(current_b,increaseDir,1),[order]])
                    
                except:
                    try:
                        am1bm1TermKinetic = np.array([hF.decrementDir(current_a,increaseDir,1),[hF.decrementDir(current_b,increaseDir,1)],[order]])
                    except:
                        am1bm1TermKinetic = np.array([[hF.decrementDir(current_a,increaseDir,1)],hF.decrementDir(current_b,increaseDir,1),[order]])    
                
                try:
                    am2TermOverlap = np.array([hF.decrementDir(current_a,increaseDir,2),current_b,[order]])
                except:
                    am2TermOverlap = np.array([[hF.decrementDir(current_a,increaseDir,2)],current_b,[order]])
                
        #-------------------------------------------------------------------------------------------        
                u = ('x','y','z')
                
                lhTerm = 'T' + hF.ind2str(lhsKinetic) + ' = '
                firstTerm = 'PA' + u[increaseDir] + '*T'  + hF.ind2str(PATermKinetic)
                    
                if current_a[increaseDir]-1 > 0:
                    secondTermKinetic = ' + ' + str(current_a[increaseDir]-1) + '*oo2p*(T' + hF.ind2str(am2TermKinetic) + ')'
                else:
                    secondTermKinetic = ''
                if current_b[increaseDir] > 0: #Not the same as for index_a
                    thirdTermKinetic = ' + ' + str(current_b[increaseDir]) + '*oo2p*(T' + hF.ind2str(am1bm1TermKinetic) + ')'
                else:
                    thirdTermKinetic = ''
    
                if current_a[increaseDir]-1 > 0: #Not the same as for index_a
                    fifthTermOverlap = ' -2*abop*' + str(current_a[increaseDir]) + '*oo2p*(S' + hF.ind2str(am2TermOverlap) + ')'
                else:
                    fifthTermOverlap = ''
                    
                fourthTerm = '+ 2*abop*S' + hF.ind2str(TermOverlap)
                stringOutput = lhTerm + firstTerm + secondTermKinetic + thirdTermKinetic + fourthTerm + fifthTermOverlap + ';\n'
                OSfile.write(stringOutput)
            elif indexNbr == 2:
                increaseDir = hF.findIncreaseDir(current_b)
        
    #            print('current_a = ' + str(current_a))        
    #            print('increasedir = ' + str(increaseDir))
    #            print(hF.decrementDir(current_a,increaseDir,1))
    #            print('current_b = ' + str(current_b))       
        
        #This part contains the actual recursion relations (i.e. the logic)--------------------------
                lhsKinetic = np.array([current_a,current_b,[order]])
                PBTermKinetic = np.array([current_a,hF.decrementDir(current_b,increaseDir,1),[order]])
                TermOverlap = np.array([current_a,current_b,[order]])
                #WQTerm = np.array([current_a,hF.decrementDir(current_b,increaseDir,1),[order+1]])
        #        print(decrementDir(current,increaseDir,2))
                try:
                    bm2TermKinetic = np.array([current_a,hF.decrementDir(current_b,increaseDir,2),[order]])
                except:
                    bm2TermKinetic = np.array([current_a,[hF.decrementDir(current_b,increaseDir,2)],[order]])
    
                try:
                    am1bm1TermKinetic = np.array([hF.decrementDir(current_a,increaseDir,1),hF.decrementDir(current_b,increaseDir,1),[order]])
                    
                except:
                    try:
                        am1bm1TermKinetic = np.array([hF.decrementDir(current_a,increaseDir,1),[hF.decrementDir(current_b,increaseDir,1)],[order]])
                    except:
                        am1bm1TermKinetic = np.array([[hF.decrementDir(current_a,increaseDir,1)],hF.decrementDir(current_b,increaseDir,1),[order]])    
                
                try:
                    bm2TermOverlap = np.array([current_a,hF.decrementDir(current_b,increaseDir,2),[order]])
                except:
                    bm2TermOverlap = np.array([current_a,[hF.decrementDir(current_b,increaseDir,2)],[order]])
                
        #-------------------------------------------------------------------------------------------        
                u = ('x','y','z')
                
                lhTerm = 'T' + hF.ind2str(lhsKinetic) + ' = '
                firstTerm = 'PB' + u[increaseDir] + '*T'  + hF.ind2str(PBTermKinetic)
                    
                if current_b[increaseDir]-1 > 0:
                    secondTermKinetic = ' + ' + str(current_b[increaseDir]-1) + '*oo2p*(T' + hF.ind2str(bm2TermKinetic) + ')'
                else:
                    secondTermKinetic = ''
                if current_a[increaseDir] > 0: #Not the same as for index_a
                    thirdTermKinetic = ' + ' + str(current_a[increaseDir]) + '*oo2p*(T' + hF.ind2str(am1bm1TermKinetic) + ')'
                else:
                    thirdTermKinetic = ''
    
                if current_b[increaseDir]-1 > 0: #Not the same as for index_a
                    fifthTermOverlap = ' -2*abop*' + str(current_b[increaseDir]-1) + '*oo2p*(S' + hF.ind2str(bm2TermOverlap) + ')'
                else:
                    fifthTermOverlap = ''
                    
                fourthTerm = '+ 2*abop*S' + hF.ind2str(TermOverlap)
                stringOutput = lhTerm + firstTerm + secondTermKinetic + thirdTermKinetic + fourthTerm + fifthTermOverlap + ';\n'
                OSfile.write(stringOutput)
            #print(stringOutput)
    #        print('lhs = '+ str(lhs))
    #        print('PATerm = '+ str(PATerm))
    #        print('WPTerm = '+ str(WPTerm))    

    return stringOutput

FileName = 'dummyTrial.m'

doOneEl(2,1,101,2,FileName)

