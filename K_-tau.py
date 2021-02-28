# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 10:24:12 2018

@author: badie
"""

from function import*
#############################   
GrotK_tau = zeros((ntheta, nk, N, 2, 2), dtype=complex) 
for l in range(0,ntheta):
    exk = cmath.exp(2*I*theta[l])   
    GrotK_tau[l] = GGrotation(Gk_tau,exk)
print GrotK_tau
###      