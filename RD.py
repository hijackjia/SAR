# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 09:11:06 2020

@author: jiash
"""
N_interp=8

import numpy as np

def interp_sinc(ys,x,Ni):
# 样本ys=f(0:N-1)
# sinc插值核长度 Ni
# 返回 y=f(x)
    N=np.size(x) #N=length(x)
    y=np.zeros(np.shape(x),dtype=complex) #y=zeros(size(x))
    ni=np.ceil(Ni/2)
    for i in range(N): # for i=1:N
        n=np.floor(x[i])
        j=np.arange(np.max([0,n-ni+1]),np.min([N-1,n+ni])+1) 
        #j=max(0,n-ni+1):min(N-1,n+ni)
        y[i]=np.dot(np.sinc(x[i]-j),ys[j.astype(int)])
    return y

def RD(sr,dt_r,dt_a,t0,RR0,VR,F0,Kr,Xmin,Xmax,Zmin,Zmax):
    c=3e8
    jpi=1j*np.pi
    #
    R_ref=np.sqrt((RR0[0]-(Xmin+Xmax)/2)**2+(RR0[2]-(Zmin+Zmax)/2)**2) # 参考距离
    lamda=c/F0
    N_r,N_a=np.shape(sr)
    tr=t0+np.arange(0,N_r)*dt_r
    #
    S=np.fft.fft2(sr)
    S=np.fft.fftshift(S)
    #
    n_r=np.ceil(N_r/2)
    fr=np.arange(n_r-N_r,n_r)/(dt_r*N_r) #fr=(n_r-N_r:n_r-1).'/(dt_r*N_r);
    n_a=np.ceil(N_a/2)
    fa=np.arange(n_a-N_a,n_a)/(dt_a*N_a) #fa=(n_a-N_a:n_a-1)/(dt_a*N_a);
    #
    D=np.sqrt(1-0.25*np.square(lamda*fa/VR))
    r_Km=0.5*np.divide(c*R_ref*np.square(fa),VR**2*F0**3*np.power(D,3)) # Z
    # r_Km=0.5*(c*R_ref*fa.^2)./(VR^2*F0^3*D.^3);
    r_Km=(1-Kr*r_Km)/Kr
    # 距离压缩
    aa=np.matmul(np.reshape(np.square(fr),(-1,1)),np.reshape(r_Km,(1,-1)))
    S=np.multiply(S,np.exp(jpi*aa))    # S=S.*exp(1i*pi*fr.^2*r_Km);
    del r_Km
    S=np.fft.ifftshift(S,axes=0)   #S=ifftshift(S,1);
    S=np.fft.ifft(S,None,axis=0)   #S=ifft(S,[],1);
    # RCMC
    for ja in range(N_a):
        tr1=(tr/D[ja]-t0)/dt_r;
        S[:,ja]=interp_sinc(S[:,ja],tr1,N_interp)
    # 方位压缩
    aa=np.matmul(np.reshape(tr,(-1,1)),np.reshape(D,(1,-1)))
    S=np.multiply(S,np.exp(2*jpi*F0*aa))#S=S.*exp(1i*2*pi*tr*F0*D);%S=S.*exp(1i*2*pi*tr*F0*(D-1));
    S=np.fft.ifftshift(S,axes=1) #S=ifftshift(S,2);
    s=np.fft.ifft(S,None,axis=1) #s=ifft(S,[],2);
    return s
