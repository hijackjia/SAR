# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 16:00:15 2020

@author: MyPC
"""
import numpy as np

def wK(sr,dt_r,dt_a,t0,RR0,VR,F0,Kr,Xmin,Xmax,Zmin,Zmax):
    c=3e8 
    jpi=1j*np.pi
    #
    R_ref=np.sqrt((RR0[0]-(Xmin+Xmax)/2)**2+(RR0[2]-(Zmin+Zmax)/2)**2)  #  参考距离
    #
    S=np.fft.fft2(sr) #S=fft2(sr) 
    N_r,N_a=np.shape(sr)
    df_r=1/(dt_r*N_r)
    df_a=1/(dt_a*N_a)
    S=np.fft.fftshift(S) #S=fftshift(S) 
    n_r=np.ceil(N_r/2)
    n0_r=n_r-N_r 
    nfr=np.arange(n0_r,n_r) #nfr=(n0_r:n_r-1).' 
    fr=nfr*df_r
    n_a=np.ceil(N_a/2) 
    n0_a=n_a-N_a 
    nfa=np.arange(n0_a,n_a) #nfa=n0_a:n_a-1
    # S=S.*(exp(-1i*2*pi*fr*t0)*ones(1,N_a)) 
    for ja in range(N_a): #ja=1:N_a
        fd2=0.25*(c*nfa[ja]*df_a/VR)**2
        a0=R_ref*np.sqrt(np.square(F0+fr)-fd2)/c
        aa=np.square(fr)/Kr
        S[:,ja]=np.multiply(S[:,ja],np.exp(jpi*(4*a0+aa)))
#        S(:,ja)=S(:,ja).*exp(1i*(4*pi*R_ref*sqrt((F0+fr).^2-fd2)/c...
#            +pi*fr.^2/Kr)) 
        fr_i=np.sqrt(np.square(F0+fr)+fd2)-F0   #  Stolt插值点
        #fr_i=sqrt((F0+fr).^2+fd2)-F0 
        if(nfa[ja]!=0):
            S[:,ja]=np.interp(fr_i,np.append(fr,fr_i[N_r-1]),np.append(S[:,ja],0))
            #S[:,ja]=interp1([fr;fr_i(N_r)],[S(:,ja);0],fr_i) 
        S[:,ja]=np.multiply(S[:,ja],np.exp(-jpi*4*R_ref*fr/c))
    s=np.fft.ifft2(np.fft.ifftshift(S)) #s=ifft2(ifftshift(S)) 
    return s