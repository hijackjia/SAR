# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 10:38:33 2020

@author: jiash
"""
from wK import wK
from CS import CS
from RD import RD
from pickle import dump
from file_s0 import load_s0
import matplotlib.pyplot as plt
import numpy as np

c=3e8
s0,dt_r,dt_a,t0,RR0,VR,F0,Kr,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax=load_s0()
#
s=wK(s0,dt_r,dt_a,t0,RR0,VR,F0,Kr,Xmin,Xmax,Zmin,Zmax)
#s=CS(s0,dt_r,dt_a,t0,RR0,VR,F0,Kr,Xmin,Xmax,Zmin,Zmax)
#s=RD(s0,dt_r,dt_a,t0,RR0,VR,F0,Kr,Xmin,Xmax,Zmin,Zmax)
#
N_r,N_a=np.shape(s)
r=c*dt_r*np.arange(0,N_r)/2 #r=c*dt_r*(0:size(s,1)-1)/2;
a=np.linspace(Ymin,Ymax,num=N_a) #y=linspace(Ymin,Ymax,size(s,2));
#
A,R=np.meshgrid(a,r)
plt.contour(A,R,np.abs(s)) #contour(r,y,abs(s)');
plt.axis('equal')  #axis equal;
plt.xlabel('azimuth');
plt.ylabel('range');
plt.show()
with open('SAR_image.pk','wb') as f: #save('SAR_image','s');
    dump(s,f)
"""
Spyder使用单独窗口画图
菜单Tools -> Preferences：
IPython Console -> Graphics，把 Backend 设置为 Qt5
重新启动 IPython Console
"""