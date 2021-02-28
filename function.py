# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 10:24:12 2018

@author: badie
"""

from variable import*


def toOmega1(ex):
    return (beta/N)*(dot(ex,Fmatrix_Bosonic))
 
def toTau1(ex):  
    return (1/beta)*(dot(inverse_Fmatrix_Bosonic,ex))

############################################
def ToTau(ex):
    list11 = zeros((N), dtype=complex)  #(N Omega)
    list12 = zeros((N), dtype=complex)
    list21 = zeros((N), dtype=complex)
    Flist11 = zeros((N), dtype=complex)
    Flist12 = zeros((N), dtype=complex)
    Flist21 = zeros((N), dtype=complex)
    F = zeros((N,2,2), dtype=complex)
    for i in range(0,N):
        list11[i] = ex[i,0 ,0]
        list12[i] = ex[i,0, 1]
        list21[i] = ex[i,1 ,0]
        list11 = array(list11)
        list12 = array(list12)
        list21 = array(list21)
       
    Flist11 = (1/beta)*dot(inverse_Fmatrix_fermionic,list11)
    Flist12 = (1/beta)*dot(inverse_Fmatrix_fermionic,list12)
    Flist21 = (1/beta)*dot(inverse_Fmatrix_fermionic,list21)
    for i in range(0,N):
        F[i,0,0] = Flist11[i]
        F[i,0,1] = Flist12[i]
        F[i,1,0] = Flist21[i]
        F[i,1,1] = Flist11[i]       
    return F
#####################
def ToOmega(ex):
    list11 = zeros((N), dtype=complex)  #(N tau)
    list12 = zeros((N), dtype=complex)
    list21 = zeros((N), dtype=complex)
    Flist11 = zeros((N), dtype=complex)
    Flist12 = zeros((N), dtype=complex)
    Flist21 = zeros((N), dtype=complex)
    F = zeros((N,2,2), dtype=complex)
    for i in range(0,N):
        list11[i] = ex[i,0 ,0]
        list12[i] = ex[i,0, 1]
        list21[i] = ex[i,1 ,0]
        list11 = array(list11)
        list12 = array(list12)
        list21 = array(list21)
       
    Flist11 = (beta/N)*dot(list11,Fmatrix_fermionic)
    Flist12 = (beta/N)*dot(list12, Fmatrix_fermionic)
    Flist21 = (beta/N)*dot(list21, Fmatrix_fermionic)
    for i in range(0,N):
        F[i,0,0] = Flist11[i]
        F[i,0,1] = Flist12[i]
        F[i,1,0] = Flist21[i]
        F[i,1,1] = Flist11[i]       
    return F
############
def TrDotProduct(tensor1, tensor2):
    product = dot(tensor1, tensor2)
    return trace(product)
##################
def coth(x):
    return (cmath.exp(2*x)+1)/(cmath.exp(2*x)-1)
###################
def G0inverse_F(Ek):
    G0 = zeros((nk, N, 2, 2), dtype=complex) #N omega and nk 
    #w1->(2*2 matrix)
    #................
    #................
    #wN->(2*2 matrix)
    for i in range(0,nk):
        for j in range(0,N):
            G0[i,j,0,0] = 1j*OmegaFermionic[j]+mu
            G0[i,j,0,1] = -Ek[i]
            G0[i,j,1,0] = -Ek[i]
            G0[i,j,1,1] = 1j*OmegaFermionic[j]+mu
    return G0
###############
def G0inverse_B(Ek):
    G0 = zeros((nk, N, 2, 2), dtype=complex) #N omega and nk 
    #w1->(2*2 matrix)
    #................
    #................
    #wN->(2*2 matrix)
    for i in range(0,nk):
        for j in range(0,N):
            G0[i,j,0,0] = 1j*OmegaBosonic[j]+mu
            G0[i,j,0,1] = -Ek[i]
            G0[i,j,1,0] = -Ek[i]
            G0[i,j,1,1] = 1j*OmegaBosonic[j]+mu
    return G0
##################
def Gtrue(Ginverse):
    G = zeros((nk,N, 2, 2), dtype=complex) #N omega   
    ev = zeros((nk,N, 2, 2), dtype=complex) #N omega
    e = zeros((nk,N, 1, 2), dtype=complex) #N omega 
    D = zeros((nk,N, 2, 2), dtype=complex) #N omega 
    for i in range(0,nk):
        for j in range(0,N):
            e[i,j], ev[i,j] = linalg.eig(eta*Ginverse[i,j])
            D[i,j,0,0] = eta*coth(e[i,j,0,0])
            D[i,j,0,1] = 0. 
            D[i,j,1,0] = 0.
            D[i,j,1,1] = eta*coth(e[i,j,0,1])       
            G[i,j] = dot(ev[i,j],dot(D[i,j],inv(ev[i,j]))) 
    return G
######
def GG(kprime, Gktau):
    Gkq = zeros((nk,N,2,2), dtype=complex)
    if kprime > max(k).real or kprime < min(k).real:
        Gkq = zeros((nk,N,2,2), dtype=complex)#N is Length TauValues
    else:
        for i in range(1,nk):
            if kprime == k[i]:
                Gkq[i] = Gktau[i]
            elif kprime == k[0]:
                Gkq[0] = Gktau[0]
            elif k[i-1] < kprime < k[i]:
                for m in range(0, N):  #N is Length TauValues
                    x = [k[i-1], k[i]]
                    y11 = [Gktau[i-1,m,0,0],Gktau[i,m,0,0]]
                    Gkq.real[i-1,m,0,0] = interp(kprime, x, y11).real
                    Gkq.imag[i-1,m,0,0] = interp(kprime, x, y11).imag

                    y12 = [Gktau[i-1,m,0,1],Gktau[i,m,0,1]]
                    Gkq.real[i-1,m,0,1] = interp(kprime, x, y12).real
                    Gkq.imag[i-1,m,0,1] = interp(kprime, x, y12).imag
                    
                    y21 = [Gktau[i-1,m,1,0],Gktau[i,m,1,0]]
                    Gkq.real[i-1,m,1,0] = interp(kprime, x, y21).real
                    Gkq.imag[i-1,m,1,0] = interp(kprime, x, y21).imag

    return Gkq
##########                
def GGrotation(G, ex):
#    G = zeros((nk,N,2,2), dtype=complex)#N is Length TauValues
    for i in range(0,nk):
        for j in range(0,N):
            G[i,j,0,0] = G[i,j,0,0]
            G[i,j,0,1] = ex*G[i,j,0,1]
            G[i,j,1,0] = conjugate(ex)*G[i,j,1,0]
            G[i,j,1,1] = G[i,j,0,0]
    return G
##################    
def Gkomega_F(sigma):
    for i in range(0,nk):
        return Gtrue(G0inverse_F(Ek)-sigma)
#################      
def Gkomega_B(sigma):
    for i in range(0,nk):
        return Gtrue(G0inverse_B(Ek)-sigma)
###############
