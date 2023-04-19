# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 14:53:54 2020

@author: MyPC
"""
import numpy as np

def CS(sr,dt_r,dt_a,t0,RR0,VR,F0,Kr,Xmin,Xmax,Zmin,Zmax):
    c=3e8 
    jpi=1j*np.pi
    #
    R_ref=np.sqrt((RR0[0]-(Xmin+Xmax)/2)**2+(RR0[2]-(Zmin+Zmax)/2)**2)  # 参考距离
    fa_ref=0  #参考方位频率
    lamda=c/F0 
    N_r,N_a=np.shape(sr)
    tr=t0+np.arange(0,N_r)*dt_r # tr=t0+(0:N_r-1).'*dt_r;
    n_r=np.ceil(N_r/2)
    fr=np.arange(n_r-N_r,n_r)/(dt_r*N_r) #fr=(n_r-N_r:n_r-1).'/(dt_r*N_r);
    n_a=np.ceil(N_a/2)
    fa=np.arange(n_a-N_a,n_a)/(dt_a*N_a) #fa=(n_a-N_a:n_a-1)/(dt_a*N_a);
    #
    S=np.fft.fft(sr,None,axis=1) #S=fft(sr,[],2) 
    S=np.fft.fftshift(S,axes=1) #S=fftshift(S,2)
    #
    D=np.sqrt(1-0.25*np.square(lamda*fa/VR)) #D=sqrt(1-0.25*(lamda*fa/VR).^2) 
    D_ref=np.sqrt(1-0.25*(lamda*fa_ref/VR)**2) #D_ref=sqrt(1-0.25*(lamda*fa_ref/VR).^2)
    Km=0.5*np.divide(c*R_ref*np.square(fa),VR**2*F0**3*np.power(D,3))
    #Km=0.5*(c*R_ref*fa.^2)./(VR^2*F0^3*D.^3) 
    Km=Kr*np.reciprocal(1-Kr*Km) #Km=Kr./(1-Kr*Km) 
    #
    a0=np.multiply(Km,D_ref*np.reciprocal(D)-1)
    a0=np.tile(a0,(N_r,1))
    aa=2*R_ref*np.reciprocal(c*D)
    aa=np.tile(np.reshape(tr,(-1,1)),N_a)-np.tile(aa,(N_r,1))
    #
    S=np.multiply(S,np.exp(jpi*np.multiply(a0,np.square(aa))))  #变标
#    S=S.*exp(1i*pi*(ones(N_r,1)*(Km.*(D_ref./D-1)))...
#        .*(tr*ones(1,N_a)-ones(N_r,1)*(2*R_ref./(c*D))).^2)  #变标
    S=np.fft.fft(S,None,axis=0)# S=fft(S,[],1) 
    S=np.fft.fftshift(S,axes=0) #S=fftshift(S,1)
    #
    a0=np.reshape(np.square(fr),(-1,1))
    a0=np.matmul(a0,np.reshape(np.divide(D,D_ref*Km),(1,-1)))
    aa=(np.reciprocal(D)-1/D_ref)/c
    aa=R_ref*np.matmul(np.reshape(fr,(-1,1)),np.reshape(aa,(1,-1)))
    #
    S=np.multiply(S,np.exp(jpi*(a0+4*aa))) 
#    S=S.*exp(1i*pi*fr.^2*(D./(D_ref*Km))... # 距离压缩+SRC
#        +1i*4*pi*R_ref*fr*(1./D-1/D_ref)/c)  # 一致RCMC
    S=np.fft.ifftshift(S,axes=0)#S=ifftshift(S,1) 
    S=np.fft.ifft(S,None,axis=0)#S=ifft(S,[],1) 
    #方位压缩
    tr=np.reshape(tr,(-1,1))
    a0=F0*np.matmul(tr,np.reshape(D,(1,-1)))
    aa=np.multiply(Km,1-D/D_ref)
    aa=np.divide(aa,np.square(D))
    aa=np.matmul(np.square(tr-2*R_ref/c),np.reshape(aa,(1,-1)))
    #
    S=np.multiply(S,np.exp(jpi*(2*a0-aa)))
    #S=S.*exp(1i*2*pi*tr*F0*(D-1)...
#    S=S.*exp(1i*2*pi*tr*F0*D...
#        -1i*pi*(tr-2*R_ref/c).^2*(Km.*(1-D/D_ref)./D.^2)) 
    S=np.fft.ifftshift(S,axes=1) #S=ifftshift(S,2) 
    s=np.fft.ifft(S,None,axis=1) #s=ifft(S,[],2) 
    return s