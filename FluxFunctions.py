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
    numNodes = len(rho)
    
    dP = np.column_stack([(gamma-1)*en, np.zeros([numNodes,2]), (gamma-1)*rho])
    
    return dP


def df12dU(rho,v1,v2,en):
    numNodes = len(rho)
    
    from scipy.sparse import csr_matrix
    
    df1rho,df1v1,df1v2,df1E,\
        df2rho,df2v1,df2v2,df2E = dfAll(rho,v1,v2,en)
    

    #number of entries in Jacobian = num vars (rho - e) squared 
    #times numNodes = 4*4*numNodes
    numEntries = 4*4*numNodes
    
    rows = np.arange(0,numEntries+1,4)
    cols = np.tile( np.arange(4), 4*numNodes) + np.tile( np.arange(numNodes), 4)
    vals1 = np.column_stack([df1rho,df1v1,df1v2,df1E]).reshape([-1])
    vals2 = np.column_stack([df2rho,df2v1,df2v2,df2E]).reshape([-1])
    
    
    df1dU = csr_matrix( (vals1,rows,cols), shape=(4*numNodes, 4*numNodes) )
    df2dU = csr_matrix( (vals2,rows,cols), shape=(4*numNodes, 4*numNodes) )
    
    return df1dU, df2dU



#This might be done now
def dfAll(rho,v1,v2,en):
    numNodes = len(rho)

    P = pressure(rho,en)
    dP = gradP(rho,en)
    
    allZero = np.zeros(numNodes)
    
    E = en + .5*( v1*v1 + v2*v2 )
    dE = np.column_stack([allZero, v1, v2, np.ones(numNodes)])
    
    
    df1_rho = np.column_stack([v1, rho, allZero, allZero])
    df1_v1 = np.column_stack([v1*v1, 2*rho*v1, allZero, allZero ]) + dP
    df1_v2 = np.column_stack([v2*v1, rho*v2, rho*v1, allZero]) + dP
    df1_E = np.column_stack([allZero, (rho*E+P), allZero, allZero]) + v1*( np.column_stack([E, allZero, allZero, allZero]) + rho*dE + dP )
    
    
    df2_rho = np.column_stack([v2, allZero, rho, allZero])
    df2_v1 = np.column_stack([v1*v2, rho*v2, rho*v1, allZero ]) + dP
    df2_v2 = np.column_stack([v2*v2, allZero, 2*rho*v2, allZero ]) + dP
    df2_E = np.column_stack([allZero, allZero, (rho*E+P), allZero]) + v2*( np.column_stack([E, allZero, allZero, allZero]) + rho*dE + dP ) 
    
    

    return df1_rho,df1_v1,df1_v2,df1_E, df2_rho,df2_v1,df2_v2,df2_E




def F12(rho,v1,v2,en):

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



