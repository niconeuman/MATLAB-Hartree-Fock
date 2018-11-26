#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  9 09:39:52 2018

@author: nneuman

This routine generates code for horizontal recursion relations

"""

import numpy as np


#This code works correctly-----------------------------------------------------
def shell(L):
    
    kmax = (L+1)*(L+2)/2
    kmax = int(kmax)
    shellInd = (np.zeros((kmax,3))).astype(int)
    
    t = 0
#    print(shellInd)
    for k in range(0,L+1):
        for m in range(0,k+1):
            shellInd[t][0] = L-k
            shellInd[t][1] = k-m
            shellInd[t][2] = m
#            print(shellInd[t])
            t = t + 1
#    print(shellInd)
    return shellInd

#------------------------------------------------------------------------------

#shellInd = shell(4)

def decrementDir(current,increaseDir,amount):
    one = np.array([[1,0,0],
                    [0,1,0],
                    [0,0,1]])
    
    result = current-amount*one[increaseDir]
#    print(result)
    if (all(result >= 0)):
        pass
    else:
        result = 0
     
    return result

def stringType(L):
    
    if (L == 0):
        strType = 's'
    elif (L == 1):
        strType = 'p'
    elif (L == 2):
        strType = 'd'
    elif (L == 3):
        strType = 'f'   
    elif (L == 4):
        strType = 'g'   
    elif (L == 5):
        strType = 'h'   
    elif (L == 6):
        strType = 'i'   
    elif (L == 7):
        strType = 'k'   
    elif (L == 8):
        strType = 'l'   
    elif (L == 9):
        strType = 'm'   
    elif (L == 10):
        strType = 'n'
    elif (L == 11):
        strType = 'o'  
    elif (L == 12):
        strType = 'q'  
    
    return strType


def ind2str(term):
    #form of term: [array([ 0.,  0.,  1.]) list([0])]
    if np.size(term) == 2:
        order = term[1][0]
        integral = term[0]
        
        L = np.sum(integral)
        
        strType = stringType(L)               
        
        u = ('x','y','z')
    #    print(integral)
        strCart = ''
        for dir in range(0,3):
    #        print(integral[dir])
            if integral[dir] > 0:
                if integral[dir] < 6:
                    for k in range(integral[dir]):
                        strCart = strCart + u[dir]
                else:
                    strCart = strCart + u[dir] + str(integral[dir])
        returnStr = strType + strCart + '_'            
    
    elif np.size(term) == 3:
        
        order = term[2][0]
        returnStr = ''
        for nshell in range(2):
            integral = term[nshell]
            
            L = np.sum(integral)
            
            strType = stringType(L)                
            
            u = ('x','y','z')
        #    print(integral)
            strCart = ''
            if np.size(integral) > 1:
                for dir in range(0,3):
            #        print(integral[dir])
                    if integral[dir] > 0:
                        if integral[dir] < 6:
                            for k in range(integral[dir]):
                                strCart = strCart + u[dir]
                        else:
                            strCart = strCart + u[dir] + str(integral[dir])
            returnStr = returnStr + strType + strCart + '_'
    elif np.size(term) == 4:
        
        order = term[3][0]
        returnStr = ''
        for nshell in range(3):
            integral = term[nshell]
            
            L = np.sum(integral)
            
            strType = stringType(L)                
            
            u = ('x','y','z')
        #    print(integral)
            strCart = ''
            if np.size(integral) > 1:
                for dir in range(0,3):
            #        print(integral[dir])
                    if integral[dir] > 0:
                        if integral[dir] < 6:
                            for k in range(integral[dir]):
                                strCart = strCart + u[dir]
                        else:
                            strCart = strCart + u[dir] + str(integral[dir])
            returnStr = returnStr + strType + strCart + '_'
            
    elif np.size(term) == 5:
        
        order = term[4][0]
        returnStr = ''
        for nshell in range(4):
            integral = term[nshell]
            
            L = np.sum(integral)
            
            strType = stringType(L)                
            
            u = ('x','y','z')
        #    print(integral)
            strCart = ''
            if np.size(integral) > 1:
                for dir in range(0,3):
            #        print(integral[dir])
                    if integral[dir] > 0:
                        if integral[dir] < 6:
                            for k in range(integral[dir]):
                                strCart = strCart + u[dir]
                        else:
                            strCart = strCart + u[dir] + str(integral[dir])
            returnStr = returnStr + strType + strCart + '_'            
#        print(returnStr)
    return returnStr + str(order)

def findIncreaseDir(current):
    if (all(current) != 0):         #If there are no zero indexes
        increaseDir = np.where(current == current.min())   #current has to be a numpy array
        increaseDir = increaseDir[0][0]                        #This had more than one dimension, with the others being empty
#            print(str(increaseDir)+ 'first')
        if np.size(increaseDir) > 1:
            increaseDir = increaseDir[-1]               #If two or more directions have the same index, increment the last one (y or z) 
#                print(increaseDir[0])
    else:
        modifiedCurrent = np.zeros(np.size(current))
#            print(modifiedCurrent)
        nind = 3
        for el in range(0,3):
            if current[el] == 0:
                modifiedCurrent[el] = 1000  #Large number which could never be an angular momentum value
                nind = nind - 1
            else:
                modifiedCurrent[el] = current[el]
        if (nind == 1):
            increaseDir = np.where(modifiedCurrent == modifiedCurrent.min())
            increaseDir = increaseDir[0][0]
        else:                 
            increaseDir = np.where(modifiedCurrent == modifiedCurrent.min())   #current has to be a numpy array
            increaseDir = increaseDir[0][0]                        #This had more than one dimension, with the others being empty
#                print(str(increaseDir)+ 'second')
            if np.size(increaseDir) > 1:
                increaseDir = increaseDir[-1]               #If two or more directions have the same index, increment the last one (y or z) 
    return increaseDir


def dohrr(La,Lc,order):
    shell_a = shell(La)
    shellInd = shell(Lc)
#    print(shellInd)
    for ta in range(np.size(shell_a,0)):
        current_a = shell_a[ta]
        for t in range(np.size(shellInd,0)):
            current = shellInd[t]
    
    
            increaseDir = findIncreaseDir(current)
    
    #        print('current = ' + str(current))        
    #        print('increasedir = ' + str(increaseDir))
    #        print(decrementDir(current,increaseDir,1))
    
    #This part contains the actual recursion relations (i.e. the logic)--------------------------
            lhs = np.array([current_a,current,[order]])
            ABTerm = np.array([current_a,decrementDir(current,increaseDir,1),[order]])
            Ap1BTerm = np.array([decrementDir(current_a,increaseDir,-1),decrementDir(current,increaseDir,1),[order]])

#            print(am1cm1Term)
                #This except case will not be executed normally, but if it were executed, 
                #it would later call ind2str with a [0] array instead of a 3-element array.
                #Should be handled
    
    #        am2Term_2 = np.array([decrementDir(current,increaseDir,2),[order+1]])
    #-------------------------------------------------------------------------------------------        
            u = ('x','y','z')
            
            lhTerm = ind2str(lhs) + ' = '
            firstTerm = 'AB' + u[increaseDir] + '*'  + ind2str(ABTerm)
            secondTerm = ind2str(Ap1BTerm)      
          
            print(lhTerm + firstTerm + ' + ' + secondTerm)
    #        print('lhs = '+ str(lhs))
    #        print('PATerm = '+ str(PATerm))
    #        print('WPTerm = '+ str(WPTerm))

    return current

def dohrr2(La,Lb,Lc,Ld,indexNbr,order):
    shellInda = shell(La)
    shellIndb = shell(Lb)
    shellIndc = shell(Lc)
    shellIndd = shell(Ld)
#    print(shellInd)
    for ta in range(np.size(shellInda,0)):
        current_a = shellInda[ta]
        for tb in range(np.size(shellIndb,0)):
            current_b = shellIndb[tb]
            for tc in range(np.size(shellIndc,0)):
                current_c = shellIndc[tc]
                for td in range(np.size(shellIndd,0)):
                    current_d = shellIndd[td]

                    if (indexNbr == 2):
                        current = current_b #just for return purposes
                        increaseDir = findIncreaseDir(current_b)
                    
                        #This part contains the actual recursion relations (i.e. the logic)--------------------------
                        lhs = np.array([current_a,current_b,current_c,current_d,[order]])
                        ABTerm = np.array([current_a,decrementDir(current_b,increaseDir,1),current_c,current_d,[order]])
                        Ap1BTerm = np.array([decrementDir(current_a,increaseDir,-1),decrementDir(current_b,increaseDir,1),current_c,current_d,[order]])
            
            #        am2Term_2 = np.array([decrementDir(current,increaseDir,2),[order+1]])
            #-------------------------------------------------------------------------------------------        
                    u = ('x','y','z')
                    
                    lhTerm = ind2str(lhs) + ' = '
                    firstTerm = 'AB' + u[increaseDir] + '*'  + ind2str(ABTerm)
                    secondTerm = ind2str(Ap1BTerm)      
                  
                    print(lhTerm + firstTerm + ' + ' + secondTerm)


    return current
    
VRR = dohrr2(2,1,1,0,2,0)
