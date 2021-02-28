# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 10:24:12 2018

@author: badie
"""

from function import*
############################# 
sigma = zeros((nk,N,2,2), dtype=complex)
Gkomega = Gkomega_F(sigma)

#print Gkomega_F

Gktau = zeros((nk,N,2,2), dtype=complex)# G(k.t)
for i in range(0,nk):
    Gktau[i] = ToTau(Gkomega[i]) 

Gk_tau = zeros((nk,N,2,2), dtype=complex)
# G(k,-t) Flensberg page(173-chapter 10)
for i in range(0, nk):
    Gk_tau[i,0] = Gktau[i,0]
    for j in range(1,N):
        Gk_tau[i,j] = Gktau[i,N-j]

############################  
KplusQ = zeros((nq, ntheta, nk),dtype=complex)
exkq = zeros((nq, ntheta, nk),dtype=complex)
GrotKplusQtau = zeros((nq,ntheta,nk,N,2,2), dtype=complex)
for j in range(0,nq):
    for l in range(0,ntheta):
        for i in range(0,nk):
            KplusQ[j,l,i] = (cmath.sqrt(k[i]**2+q[j]**2+2*k[i]*q[j]*\
            cos(theta[l]))).real
            exkq[j,l,i] = (k[i]*cmath.exp(1j*theta[l])+q[j])/\
            (k[i]*cmath.exp(-1j*theta[l])+q[j])

for j in range(0,nq):
    for l in range(0,ntheta):
        for i in range(0,nk):
            if KplusQ[j,l,i] < max(k).real or KplusQ[j,l,i] == max(k).real:
                kprime = KplusQ[j,l,i]
                ex = exkq[j,l,i]
                GrotKplusQtau[j,l] = GGrotation(GG(kprime,Gktau),ex)

#print GrotKplusQtau#[:,:,:,2]
#######
#############################   
GrotK_tau = zeros((ntheta, nk, N, 2, 2), dtype=complex) 
for l in range(0,ntheta):
    exk = cmath.exp(2*1j*theta[l])   
    GrotK_tau[l] = GGrotation(Gk_tau,exk)# in here  GGrotation==Gk_tau
#print GrotK_tau
###   
GrotKplusQtau[:,:,:,0] = -0.5*identity(2)
GrotK_tau[:,:,0] = 0.5*identity(2)
#print "GrotKplusQtau =", GrotKplusQtau
#print "GrotK_tau =", GrotK_tau
###########
Matrixproduct = zeros((nq,ntheta,nk,N,1), dtype=complex)
for i in range(0,nq):
    for j in range(0,ntheta):
        for m in range(0,N):
            for n in range(1, nk-1):  # N is numbr of tau
                Matrixproduct[i,j,n,m] = ((k[n+1]-k[n-1])*k[n]/2)\
                *trace(dot(GrotK_tau[j,n,m], GrotKplusQtau[i,j,n,m]))
#print Matrixproduct#[60,0,0,3]
###
PPi = zeros((nq,N), dtype=complex)
for i in range(0,nq):
    for m in range(0,N):
        PPi.real[i,m] = (sum(sum(Matrixproduct[i,j,n,m].real\
        for j in xrange(0,ntheta)) for n in xrange(0,nk)))
        PPi.imag[i,m] = (sum(sum(Matrixproduct[i,j,n,m].imag\
        for j in xrange(0,ntheta)) for n in xrange(0,nk)))
PPi_Qtau = (-4/(2*pi*ntheta))*PPi
#print PPi_Qtau
PPi_Qw = zeros((nq,N), dtype=complex)
for i in range(0,nq):
    PPi_Qw[i] = toOmega1(PPi_Qtau[i])
#print PPi_Qw
##
Veff_QiW = zeros((nq,N),dtype=complex)
for i in range(0,nq):
    for j in range(0,N):
        Veff_QiW[i,j] = vq[i]/(1+vq[i]*PPi_Qw[i,j])
#print Veff_QiW
##########
epsilon = zeros((nq,N),dtype=complex)
for i in range(0,nq):
    for j in range(0,N):
        epsilon[i,j] = 1+vq[i]*PPi_Qw[i,j]
        
#epsilon = (-epsilon**(-1)).imag
#x = arange(0,4)
#y= arange(0,4)
#X,Y = meshgrid(x,y)
#plt.contourf(X,Y,epsilon)
#plt.contourf(epsilon)
#plt.show()
#########
Veff_Qtau = zeros((nq,N),dtype=complex)
for i in range(0,nq):
    Veff_Qtau[i] = toTau1(Veff_QiW[i])
#print Veff_Qtau[2,3]
######
KminusQ = zeros((nq, ntheta, nk),dtype=complex)
exk_q = zeros((nq, ntheta, nk),dtype=complex)
GrotKminusQtau = ones((nq,ntheta,nk,N,2,2), dtype=complex)
for j in range(0,nq):
    for l in range(0,ntheta):
        for i in range(0,nk):
            KminusQ[j,l,i] = (cmath.sqrt(k[i]**2+q[j]**2-2*k[i]*q[j]*\
            cos(theta[l]))).real
            exk_q[j,l,i] = (-q[j]*cmath.exp(1j*theta[l])+k[i])/\
            (-q[j]*cmath.exp(-1j*theta[l])+k[i])

for j in range(0,nq):
    for l in range(0,ntheta):
        for i in range(0,nk):
            if KminusQ.real[j,l,i] < max(k).real \
            or KplusQ.real[j,l,i] == max(k).real:
                kprime = KminusQ[j,l,i]
                ex = exk_q[j,l,i]
                GrotKminusQtau[j,l] = GGrotation(GG(kprime,Gktau),ex)
#print GrotKminusQtau
#print GrotKminusQtau[:,:,:,2]
#######
SS = zeros((nq,ntheta,nk,N,2,2), dtype=complex)

for i in range(1,nq-1):
    for m in range(0,ntheta):
        for l in range(0,nk):
            for n in range(0,N):
                SS[i,m,l,n] = ((q[i+1]-q[i-1])*q[i]/2)*\
                Veff_Qtau[i,n]*GrotKminusQtau[i,m,l,n]
#print SS
sigmaktau = zeros((nk,N,2,2),dtype=complex)       
sigmakomega =  zeros((nk,N,2,2),dtype=complex)
for l in range(0,nk):
    for n in range(0,N):
        sigma.real[l,n] = (sum(sum(SS[i,m,l,n].real\
        for m in xrange(0,ntheta)) for i in xrange(0,nq)))
        sigma.imag[l,n] = (sum(sum(SS[i,m,l,n].imag\
        for m in xrange(0,ntheta)) for i in xrange(0,nq)))
sigmaktau = (-1/(2*pi*ntheta))*sigma 
for i in range(0,nk):
    sigmakomega[i] = ToOmega(sigmaktau[i])
#print sigmakomega
#print sigmaktau
#    ####################################
while True:
    Gkomega_new =  zeros((nk,N,2,2),dtype=complex)  
    Gkomega_new = Gkomega_F(sigmakomega)
    Gktau_new = zeros((nk,N,2,2), dtype=complex)# G(k.t)
    for i in range(0,nk):
        Gktau_new[i] = ToTau(Gkomega_new[i]) 
    Gk_tau_new = zeros((nk,N,2,2), dtype=complex)
    # G(k,-t) Flensberg page(173-chapter 10)
    for i in range(0, nk):
        Gk_tau_new[i,0] = Gktau_new[i,0]
        for j in range(1,N):
            Gk_tau_new[i,j] = Gktau_new[i,N-j]
##
    Gkomega_new_inverse = zeros((nk,N,2,2),dtype=complex)
    for i in range(0, nk):
        for j in range(0,N):
            Gkomega_new_inverse[i,j] = inv(Gkomega_F(sigmakomega)[i,j])    
#print Gkomega_new_inverse
###########################################################
    for j in range(0,nq):
        for l in range(0,ntheta):
            for i in range(0,nk):
                if KplusQ[j,l,i] < max(k).real or\
                KplusQ[j,l,i] == max(k).real:
                    kprime = KplusQ[j,l,i]
                    ex = exkq[j,l,i]
                    GrotKplusQtau[j,l] = GGrotation(GG(kprime,Gktau_new),ex)            
#####################################################
    for l in range(0,ntheta):
        exk = cmath.exp(2*1j*theta[l])   
        GrotK_tau[l] = GGrotation(Gk_tau_new,exk)
    GrotKplusQtau[:,:,:,0] = -0.5*identity(2)
    GrotK_tau[:,:,0] = 0.5*identity(2)
#######
    Matrixproduct = zeros((nq,ntheta,nk,N,1), dtype=complex)
    for i in range(0,nq):
        for j in range(0,ntheta):
            for m in range(0,N):
                for n in range(1, nk-1):  # N is numbr of tau
                    Matrixproduct[i,j,n,m] = ((k[n+1]-k[n-1])*k[n]/2)\
                    *trace(dot(GrotK_tau[j,n,m], GrotKplusQtau[i,j,n,m]))

    PPi = zeros((nq,N), dtype=complex)
    for i in range(0,nq):
        for m in range(0,N):
            PPi.real[i,m] = (sum(sum(Matrixproduct[i,j,n,m].real\
            for j in xrange(0,ntheta)) for n in xrange(0,nk)))
            PPi.imag[i,m] = (sum(sum(Matrixproduct[i,j,n,m].imag\
            for j in xrange(0,ntheta)) for n in xrange(0,nk)))
    PPi_Qtau = (-4/(2*pi*ntheta))*PPi
#print PPi_Qtau
    PPi_Qw = zeros((nq,N), dtype=complex)
    for i in range(0,nq):
        PPi_Qw[i] = toOmega1(PPi_Qtau[i])

##
    Veff_QiW_new = zeros((nq,N),dtype=complex)
    for i in range(0,nq):
        for j in range(0,N):
            Veff_QiW_new[i,j] = vq[i]/(1+vq[i]*PPi_Qw[i,j])

    Veff_Qtau_new = zeros((nq,N),dtype=complex)
    for i in range(0,nq):
        Veff_Qtau_new[i] = toTau1(Veff_QiW_new[i])
    Delta_Veff = zeros((nq,N),dtype=complex)   
    Delta_Veff = abs(Veff_Qtau_new - Veff_Qtau)
    
    if Delta_Veff.max().real > 10**(-30):
        Veff_Qtau = Veff_Qtau_new
        print "Delta_Veff.max=",Delta_Veff.max()
        for j in range(0,nq):
            for l in range(0,ntheta):
                for i in range(0,nk):
                    if KminusQ.real[j,l,i] < max(k).real \
                    or KplusQ.real[j,l,i] == max(k).real:
                        kprime = KminusQ[j,l,i]
                        ex = exk_q[j,l,i]
                        GrotKminusQtau[j,l] =\
                        GGrotation(GG(kprime,Gktau_new),ex)
                        
        for i in range(1,nq-1):
            for m in range(0,ntheta):
                for l in range(0,nk):
                    for n in range(0,N):
                        SS[i,m,l,n] = ((q[i+1]-q[i-1])*q[i]/2)*\
                        Veff_Qtau_new[i,n]*GrotKminusQtau[i,m,l,n]

        sigmaktau = zeros((nk,N,2,2),dtype=complex)       
        sigmakomega =  zeros((nk,N,2,2),dtype=complex)
        for l in range(0,nk):
            for n in range(0,N):
                sigma.real[l,n] = (sum(sum(SS[i,m,l,n].real\
                for m in xrange(0,ntheta)) for i in xrange(0,nq)))
                sigma.imag[l,n] = (sum(sum(SS[i,m,l,n].imag\
                for m in xrange(0,ntheta)) for i in xrange(0,nq)))
        sigmaktau = (-1/(2*pi*ntheta))*sigma 
        for i in range(0,nk):
            sigmakomega[i] = ToOmega(sigmaktau[i])
            
    else:
        break

print Gkomega_new                  
##########
#epsilon = zeros((nq,N),dtype=complex)
#for i in range(0,nq):
#    for j in range(0,N):
#        epsilon[i,j] = 1+vq[i]*PPi_Qw[i,j]