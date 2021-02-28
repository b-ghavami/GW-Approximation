# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 10:24:12 2018

@author: badie
"""

from Lib import*
#####################3
N=int(input("please enter number of frequency : "))
#N = 4
beta = 10.0
mu = 1.0
rs = 7.0
k1 = 0.0
k2 = 4.0
nk = 41
nk1 = 19
nk2 = 19
ntheta = 10
eta = beta/(2*N)
#################
lambda1 = abs(log(1-k1)-4 * log(10))
lambda2 = abs(log(k2-1) + 4 * log(10))
#
k = zeros((nk), dtype=complex)
for s in range(0, nk1+1, 1):
    k[s] = 1 - cmath.exp(-(s) * lambda1 / nk1)
    k[20] = 1.0
for s in range(0, nk2+1, 1):
    k[s + 21] = 1 + cmath.exp(-lambda1 + (s) * lambda2 / nk2)
k = array(k) # k is a matrix with one column
#print k
Ek = zeros((nk),dtype=complex)
Ek=k**2
#plt.plot(k,Ek)
#print Ek[5]
nqq = 81
qq1 = 1.0 / nqq
qq2 = 4
nq = nqq + 1
lambdaq = abs(qq2 - qq1)
qq = zeros((nq), dtype=complex)
vq = zeros((nq), dtype=complex)
for s in range(0, nq, 1):
    qq[s] = qq1 + (s) * lambdaq / nqq
    vq[s] = pi * rs / qq[s]
#    print qq
q = array(qq)
vq = array(vq)
#print q
#print vq
#plt.plot(qq,vq)
OmegaFermionic = zeros((N), dtype=complex)
for i in range(0, N):
    OmegaFermionic[i] = (2 * (i) + 1) * pi / beta
#OmegaFermionic = mat([OmegaFermionic]).T
OmegaFermionic = array(OmegaFermionic)
#print OmegaFermionic, "===="
#
OmegaBosonic = zeros((N), dtype=complex)
for i in range(0, N):
    OmegaBosonic[i] = (2 * (i) ) * pi / beta
OmegaBosonic = array(OmegaBosonic)
#print OmegaBosonic,"====" 
OmegaFB = zeros((N,N), dtype=complex)
for i in range(0,N):
    for j in range(0,N):
        OmegaFB[i,j] = OmegaFermionic[i]+OmegaBosonic[j]#WF+WB
#print OmegaFB
#
theta = zeros((ntheta), dtype=complex)
for s in range(0, ntheta):
    theta[s] = 2 * pi * (s) / ntheta + pi / ntheta
theta = array(theta)
#print theta[3,0]
#
tau = zeros((N), dtype=complex)#is a matrix with one row
for i in range(1, N):
    tau[i] = ( i * beta) / N
tau = array(tau)
#print tau

#W = array([omega])
#print TensorProduct(W,W.conj().T)
#print W.conj().T
Fmatrix_fermionic = zeros((N, N), dtype=complex)
inverse_Fmatrix_fermionic = zeros((N, N), dtype=complex)
for i in range(0,N):
    for j in range(0,N):
        Fmatrix_fermionic[i,j] = cmath.exp(tau[i]*OmegaFermionic[j]*1j)
        inverse_Fmatrix_fermionic[i,j] = cmath.exp(-tau[i]\
        *OmegaFermionic[j]*1j)
Fmatrix_fermionic = array(Fmatrix_fermionic)
inverse_Fmatrix_fermionic = array(inverse_Fmatrix_fermionic)
#inverse_Fmatrix_fermionic = mat(inv(Fmatrix_fermionic))
#inverse_Fmatrix_fermionic = inv(Fmatrix_fermionic)
#print Fmatrix_fermionic, '===', inverse_Fmatrix_fermionic
#print inverse_Fmatrix_fermionic*Fmatrix_fermionic
#######
Fmatrix_Bosonic = zeros((N, N), dtype=complex)
inverse_Fmatrix_Bosonic = zeros((N, N), dtype=complex)

for i in range(0,N):
    for j in range(0,N):
        Fmatrix_Bosonic[i,j] = cmath.exp(tau[i]*OmegaBosonic[j]*1j)
        inverse_Fmatrix_Bosonic[i,j] = cmath.exp(-tau[i]\
        *OmegaBosonic[j]*1j)
Fmatrix_Bosonic = array(Fmatrix_Bosonic)
inverse_Fmatrix_Bosonic = array(inverse_Fmatrix_Bosonic)
#print Fmatrix_Bosonic
#print "-------------------------"
#print inverse_Fmatrix_Bosonic
#print "////////////"
#print dot(Fmatrix_Bosonic,inverse_Fmatrix_Bosonic)
#--------------------------------------------------------------

