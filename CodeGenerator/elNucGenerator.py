#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 11:27:41 2018

@author: nneuman

This routine generates code for one-electron electron-Nuclear integrals
eN recursion relation
[a+1i|V|b](m) = PAi*[a|V|b](m)-PCi*[a|V|b](m+1)
                +ai*(1/2p)*([a-1i|V|b](m)-[a-1i|V|b](m+1))
                +bi*(1/2p)*([a|V|b-1i](m)-[a|V|b-1i](m+1))
"""

import numpy as np
import helperFunctions as hF


def doElNuc(La,Lb,order,indexNbr,FileName):
    OSfile = open(FileName,'a+')
    #I require something that reverts to vrr if Lc == 0
    shell_a = hF.shell(La)
    shell_b = hF.shell(Lb)
#    print(shell_a)
    
    #Electron Nuclear
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
                lhs = np.array([current_a,current_b,[order]])
                PATerm = np.array([hF.decrementDir(current_a,increaseDir,1),current_b,[order]])
                PCTerm = np.array([hF.decrementDir(current_a,increaseDir,1),current_b,[order+1]])
                
                #WQTerm = np.array([current_a,hF.decrementDir(current_b,increaseDir,1),[order+1]])
        #        print(decrementDir(current,increaseDir,2))
                try:
                    am2Term_1 = np.array([hF.decrementDir(current_a,increaseDir,2),current_b,[order]])
                    am2Term_2 = np.array([hF.decrementDir(current_a,increaseDir,2),current_b,[order+1]])
                except:
                    am2Term_1 = np.array([[hF.decrementDir(current_a,increaseDir,2)],current_b,[order]])
                    am2Term_2 = np.array([[hF.decrementDir(current_a,increaseDir,2)],current_b,[order+1]])
                try:
                    am1bm1Term_1 = np.array([hF.decrementDir(current_a,increaseDir,1),hF.decrementDir(current_b,increaseDir,1),[order]])
                    am1bm1Term_2 = np.array([hF.decrementDir(current_a,increaseDir,1),hF.decrementDir(current_b,increaseDir,1),[order+1]])
                except:
                    try:
                        am1bm1Term_1 = np.array([hF.decrementDir(current_a,increaseDir,1),[hF.decrementDir(current_b,increaseDir,1)],[order]])
                        am1bm1Term_2 = np.array([hF.decrementDir(current_a,increaseDir,1),[hF.decrementDir(current_b,increaseDir,1)],[order+1]])
                    except:
                        am1bm1Term_1 = np.array([[hF.decrementDir(current_a,increaseDir,1)],hF.decrementDir(current_b,increaseDir,1),[order]]) 
                        am1bm1Term_2 = np.array([[hF.decrementDir(current_a,increaseDir,1)],hF.decrementDir(current_b,increaseDir,1),[order+1]]) 
                
                
        #-------------------------------------------------------------------------------------------        
                u = ('x','y','z')
                
                lhTerm = 'V' + hF.ind2str(lhs) + ' = '
                firstTerm = 'PA' + u[increaseDir] + '*V'  + hF.ind2str(PATerm)
                secondTerm = ' + ' + 'PC' + u[increaseDir] + '*V'  + hF.ind2str(PCTerm)  
                if current_a[increaseDir]-1 > 0:
                    thirdTerm = ' + ' + str(current_a[increaseDir]-1) + '*oo2p*(V' + hF.ind2str(am2Term_1) + ' - V' + hF.ind2str(am2Term_2) + ')'
                else:
                    thirdTerm = ''
                if current_b[increaseDir] > 0: #Not the same as for index_a
                    fourthTerm = ' + ' + str(current_b[increaseDir]) + '*oo2p*(V' + hF.ind2str(am1bm1Term_1) + ' - V' + hF.ind2str(am1bm1Term_2) + ')'
                else:
                    fourthTerm = ''
    
                stringOutput = lhTerm + firstTerm + secondTerm + thirdTerm + fourthTerm + ';\n'
                OSfile.write(stringOutput)
                
            elif indexNbr == 2:
                increaseDir = hF.findIncreaseDir(current_b)
        
    #            print('current_a = ' + str(current_a))        
    #            print('increasedir = ' + str(increaseDir))
    #            print(hF.decrementDir(current_a,increaseDir,1))
    #            print('current_b = ' + str(current_b))       
        
        #This part contains the actual recursion relations (i.e. the logic)--------------------------
                lhs = np.array([current_a,current_b,[order]])
                PBTerm = np.array([current_a,hF.decrementDir(current_b,increaseDir,1),[order]])
                PCTerm = np.array([current_a,hF.decrementDir(current_b,increaseDir,1),[order+1]])
                
                #WQTerm = np.array([current_a,hF.decrementDir(current_b,increaseDir,1),[order+1]])
        #        print(decrementDir(current,increaseDir,2))
                try:
                    bm2Term_1 = np.array([current_a,hF.decrementDir(current_b,increaseDir,2),[order]])
                    bm2Term_2 = np.array([current_a,hF.decrementDir(current_b,increaseDir,2),[order+1]])
                except:
                    bm2Term_1 = np.array([current_a,[hF.decrementDir(current_b,increaseDir,2)],[order]])
                    bm2Term_2 = np.array([current_a,[hF.decrementDir(current_b,increaseDir,2)],[order+1]])
                try:
                    am1bm1Term_1 = np.array([hF.decrementDir(current_a,increaseDir,1),hF.decrementDir(current_b,increaseDir,1),[order]])
                    am1bm1Term_2 = np.array([hF.decrementDir(current_a,increaseDir,1),hF.decrementDir(current_b,increaseDir,1),[order+1]])
                except:
                    try:
                        am1bm1Term_1 = np.array([hF.decrementDir(current_a,increaseDir,1),[hF.decrementDir(current_b,increaseDir,1)],[order]])
                        am1bm1Term_2 = np.array([hF.decrementDir(current_a,increaseDir,1),[hF.decrementDir(current_b,increaseDir,1)],[order+1]])
                    except:
                        am1bm1Term_1 = np.array([[hF.decrementDir(current_a,increaseDir,1)],hF.decrementDir(current_b,increaseDir,1),[order]]) 
                        am1bm1Term_2 = np.array([[hF.decrementDir(current_a,increaseDir,1)],hF.decrementDir(current_b,increaseDir,1),[order+1]]) 
                
                
        #-------------------------------------------------------------------------------------------        
                u = ('x','y','z')
                
                lhTerm = 'V' + hF.ind2str(lhs) + ' = '
                firstTerm = 'PB' + u[increaseDir] + '*V'  + hF.ind2str(PBTerm)
                secondTerm = ' + ' + 'PC' + u[increaseDir] + '*V'  + hF.ind2str(PCTerm)  
                if current_b[increaseDir]-1 > 0:
                    thirdTerm = ' + ' + str(current_b[increaseDir]-1) + '*oo2p*(V' + hF.ind2str(bm2Term_1) + ' - V' + hF.ind2str(bm2Term_2) + ')'
                else:
                    thirdTerm = ''
                if current_a[increaseDir] > 0: #Not the same as for index_a
                    fourthTerm = ' + ' + str(current_a[increaseDir]) + '*oo2p*(V' + hF.ind2str(am1bm1Term_1) + ' - V' + hF.ind2str(am1bm1Term_2) + ')'
                else:
                    fourthTerm = ''
    
                stringOutput = lhTerm + firstTerm + secondTerm + thirdTerm + fourthTerm + ';\n'
                OSfile.write(stringOutput)
            #print(stringOutput)
    #        print('lhs = '+ str(lhs))
    #        print('PATerm = '+ str(PATerm))
    #        print('WPTerm = '+ str(WPTerm))    

    return stringOutput

FileName = 'dummyTrial2.m'

doElNuc(2,1,0,2,FileName)