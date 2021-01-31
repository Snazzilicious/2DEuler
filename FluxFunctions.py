#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 15:48:44 2019

@author: ian
"""

import numpy as np

gamma = 1.4


def pressure(rho,en):
    return (gamma-1)*rho*en


def gradP(rho,en):
    numVerts = len(rho)
    
    dP = np.array([(gamma-1)*en, np.zeros([numVerts,2]), (gamma-1)*rho])
    
    return dP


def df123dU(rho,v1,v2,en):

    from scipy.sparse import coo_matrix
    
    df1rho,df1v1,df1v2,df1E,\
        df2rho,df2v1,df2v2,df2E = dfAll(rho,v1,v2,en)
    

    #number of entries in Jacobian = num vars (rho - e) squared 
    #times numCells = 5*5*numCells
    numEntries = 5*5*numCells
    
    rows = np.zeros([numEntries,1])
    cols = np.zeros([numEntries,1])
    vals1 = np.zeros([numEntries,1])
    vals2 = np.zeros([numEntries,1])
    vals3 = np.zeros([numEntries,1])
    
    oneToNumCells = np.arange(numCells).transpose()
    
    #going left to right top to bottom
    for i in range(4):
        rows[oneToNumCells + i*numCells] = oneToNumCells
        cols[oneToNumCells + i*numCells] = oneToNumCells + i*numCells
        vals1[oneToNumCells + i*numCells] = df1rho[:,i]
        vals2[oneToNumCells + i*numCells] = df2rho[:,i]
        vals3[oneToNumCells + i*numCells] = df3rho[:,i]

    
    for i in range(4):
        rows[oneToNumCells + [i+5]*numCells] = oneToNumCells + numCells
        cols[oneToNumCells + [i+5]*numCells] = oneToNumCells + i*numCells
        vals1[oneToNumCells + [i+5]*numCells] = df1v1[:,i]
        vals2[oneToNumCells + [i+5]*numCells] = df2v1[:,i]
        vals3[oneToNumCells + [i+5]*numCells] = df3v1[:,i]

    
    for i in range(4):
        rows[oneToNumCells + [i+10]*numCells] = oneToNumCells + 2*numCells
        cols[oneToNumCells + [i+10]*numCells] = oneToNumCells + i*numCells
        vals1[oneToNumCells + [i+10]*numCells] = df1v2[:,i]
        vals2[oneToNumCells + [i+10]*numCells] = df2v2[:,i]
        vals3[oneToNumCells + [i+10]*numCells] = df3v2[:,i]

    
    for i in range(4):
        rows[oneToNumCells + [i+15]*numCells] = oneToNumCells + 3*numCells
        cols[oneToNumCells + [i+15]*numCells] = oneToNumCells + i*numCells
        vals1[oneToNumCells + [i+15]*numCells] = df1V3[:,i]
        vals2[oneToNumCells + [i+15]*numCells] = df2V3[:,i]
        vals3[oneToNumCells + [i+15]*numCells] = df3V3[:,i]

    
    for i in range(4):
        rows[oneToNumCells + [i+20]*numCells] = oneToNumCells + 4*numCells
        cols[oneToNumCells + [i+20]*numCells] = oneToNumCells + i*numCells
        vals1[oneToNumCells + [i+20]*numCells] = df1E[:,i]
        vals2[oneToNumCells + [i+20]*numCells] = df2E[:,i]
        vals3[oneToNumCells + [i+20]*numCells] = df3E[:,i]

    
    df1dU = coo_matrix(vals1, (rows,cols), shape=(5*numCells, 5*numCells) )
    df2dU = coo_matrix(vals2, (rows,cols), shape=(5*numCells, 5*numCells) )
    
    return df1dU, df2dU



#This is not done
def dfAll(rho,v1,v2,V3,en):
    global numCells

    P = pressure(rho,en)
    dP = gradP(rho,en)
    
    allZero = np.zeros(numCells)
    
    E = en + .5*( v1*v1 + v2*v2 )
    dE = np.array([allZero, allZero, allZero, V3, np.ones(numCells)]) + .5*GR.gradCrossFlowMag2(v1,v2,gg) 
    
    
    df1_rho = np.array([v1, rho, allZero, allZero, allZero])
    df1_v1 = np.array([v1*v1, 2*rho*v1, allZero, allZero, allZero ]) + gInv[1,1,:]*dP
    df1_v2 = np.array([v2*v1, rho*v2, rho*v1, allZero, allZero]) + gInv[2,1,:]*dP
    df1_E = np.array([allZero, (rho*E+P), allZero, allZero, allZero]) + v1*( [E, allZero, allZero, allZero, allZero] + rho*dE + dP )
    
    
    df2_rho = np.array([v2, allZero, rho, allZero, allZero])
    df2_v1 = np.array([v1*v2, rho*v2, rho*v1, allZero, allZero ]) + gInv[1,2,:]*dP
    df2_v2 = np.array([v2*v2, allZero, 2*rho*v2, allZero, allZero ]) + gInv[2,2,:]*dP
    df2_E = np.array([allZero, allZero, (rho*E+P), allZero, allZero]) + v2*( [E, allZero, allZero, allZero, allZero] + rho*dE + dP ) 
    
    

    return df1_rho,df1_v1,df1_v2,df1_E, df2_rho,df2_v1,df2_v2,df2_E




def F123(rho,v1,v2,en):

    f1_rho,f1_v1,f1_v2,f1_E,\
        f2_rho,f2_v1,f2_v2,f2_E = fAll(rho,v1,v2,en)
    
#    F1 = np.array([f1_rho, f1_v1, f1_v2, f1_V3, f1_E]).reshape(-1,1)
#    F2 = np.array([f2_rho, f2_v1, f2_v2, f2_V3, f2_E]).reshape(-1,1)
#    F3 = np.array([f3_rho, f3_v1, f3_v2, f3_V3, f3_E]).reshape(-1,1)
    
    F1 = np.concatenate([f1_rho, f1_v1, f1_v2, f1_E])
    F2 = np.concatenate([f2_rho, f2_v1, f2_v2, f2_E])
    
    return F1, F2



def fAll(rho,v1,v2,en):
    
    P = pressure(rho,en)
    
    E = en + .5*( v1*v1 + v2*v2 )

    f1_rho = rho*v1
    f1_v1 = rho*v1*v1 + P
    f1_v2 = rho*v2*v1 + P
    f1_E = (rho*E + P)*v1 
    
    f2_rho = rho*v2
    f2_v1 = rho*v1*v2 + P
    f2_v2 = rho*v2*v2 + P
    f2_E = (rho*E + P)*v2 
    
    
    return f1_rho,f1_v1,f1_v2,f1_E, f2_rho,f2_v1,f2_v2,f2_E



