# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 10:24:12 2018

@author: badie
"""

from function import*
#############################   
#KplusQ = zeros((nk+nq+ntheta),dtype=complex)
#ex = zeros((nk+nq+ntheta),dtype=complex)
GrotKplusQtau = zeros((ntheta,nq,nk,N,2,2), dtype=complex)
for i in range(0,nk):
    for j in range(0,nq):
        for l in range(0,ntheta):
            KplusQ = (cmath.sqrt(k[i]**2+q[j]**2+2*k[i]*q[j]*\
            cos(theta[l]))).real
            #
            exkq = (k[i]*cmath.exp(I*theta[l])+q[j])/\
            (k[i]*cmath.exp(-I*theta[l])+q[j])
#            print KplusQ
            if KplusQ < max(k).real or KplusQ == max(k).real:
                GG(KplusQ, Gktau)
                print GG(KplusQ, Gktau)
#            
                for m in range(0,nq):
                    GrotKplusQtau[l,m] = GGrotation(GG(KplusQ,Gktau),exkq)
#print GrotKplusQtau
########
