# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 11:20:36 2020

@author: jiash
"""
from pickle import load
from pickle import dump

def save_s0(s0,dt_r,dt_a,t0,RR0,VR,F0,Kr,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax):
    with open('SAR_s0.pk','wb') as f:
        dump(s0,f)
        dump(dt_r,f)
        dump(dt_a,f)
        dump(t0,f)
        dump(RR0,f)
        dump(VR,f)
        dump(F0,f)
        dump(Kr,f)
        dump(Xmin,f)
        dump(Xmax,f)
        dump(Ymin,f)
        dump(Ymax,f)        
        dump(Zmin,f)
        dump(Zmax,f)

def load_s0():
    with open('SAR_s0.pk','rb') as f:
        s0  =load(f)
        dt_r=load(f)
        dt_a=load(f)
        t0  =load(f)
        RR0 =load(f)
        VR  =load(f)
        F0  =load(f)
        Kr  =load(f)
        Xmin=load(f)
        Xmax=load(f)
        Ymin=load(f)
        Ymax=load(f)        
        Zmin=load(f)
        Zmax=load(f)
    return s0,dt_r,dt_a,t0,RR0,VR,F0,Kr,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax