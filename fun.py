# -*- coding: utf-8 -*-
"""
Created on Sat Jun 13 08:54:43 2020

@author: MyPC
"""
import numpy as np

def antenna_pattern(theta,theta_bw):
    return np.square(np.sinc(theta/theta_bw))

def Env_pul(t,T):
    x=np.abs(t)<=T/2
    return x.astype(float) 

def waveform(R,NT,N_r,d_w,tr,Tr,Kr,F0,Theta_bw): # R(NT,3)
    c=3e8
    jpi=1j*np.pi
    #
    d=np.linalg.norm(R,axis=1)
    dd=np.tile(np.reshape(d,(-1,1)),3)
    R=np.divide(R,dd)
    theta=np.arccos(np.matmul(R,d_w))
    #
    tr=np.tile(tr,(NT,1))
    td=2*d/c
    td=np.tile(np.reshape(td,(-1,1)),N_r)
    t=tr-td
    #
    s=Env_pul(t,Tr)
    ap=antenna_pattern(theta,Theta_bw)
    s=np.multiply(s,np.tile(np.reshape(ap,(-1,1)),N_r))
    s=np.multiply(s,np.exp(jpi*(Kr*np.square(t)-2*F0*td)))
    return s


    
