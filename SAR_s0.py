# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 13:19:22 2020

@author: MyPC
"""
Xmin=0
Xmax=5
Ymin=0
Ymax=4
Zmin=-1
Zmax=2  # 目标坐标范围(m)
#####################################################################
Dr  = 0.1      # 距离向分辨率 =c/(2*BW)
#
Da  = 0.1       # 方位向分辨率 =VR/BW_a =lamda_min/(2*sin(Theta_bw))
# 方位向带宽为最大多普勒频移 BW_a = 2*F_max*VR*sin(Theta_bw)/c
#
angle_bw = 0       # 波束中心方向与水平面（xy面）夹角(角度deg)
#     \   /|\  z
#     _\|  |
#       \  |
#        \ |
# angle_bw\|____________\  x
#                       /
####################################################################
VR  =  3       # 雷达速度 (沿正y方向)(m/s)
A_os_r = 1.2      # 距离向过采样因子(>=1)
A_os_a = 1.2      # 方位向过采样因子(>=1)
#
import numpy as np
from pickle import load
from file_s0 import save_s0
from fun import waveform
from sys import exit
#
angle_bw=angle_bw*np.pi/180
d_w =np.array([np.cos(angle_bw),0,-np.sin(angle_bw)]) # 波束中心方向矢量
c=3e8 
#
TBP = 100       # 时间带宽积(>=100)
BW  = c/(2*Dr)      # 带宽
Tr  = TBP/BW    #  脉冲时宽
Kr  = BW/Tr         #调频率
dt_r= 1/(BW*A_os_r) # 距离向采样间隔
F0  = 20*BW       # 中心频率
#
lamda=c/F0 
lamda_min=c/(F0+BW/2) 
Theta_bw = np.arcsin(0.5*lamda_min/Da) # 天线波束半角宽
BW_a=2*VR*np.sin(Theta_bw)/lamda_min  # 方位向带宽
dt_a=1/(BW_a*A_os_a)  #方位向时间步长
###################################
Rmin=np.array([Xmin,Ymin,Zmin]) 
Rmax=np.array([Xmax,Ymax,Zmax])
#D=sqrt(2)*HR*tan(Theta_bw) 
D=0.5*np.sqrt((Xmax-Xmin)**2+(Zmax-Zmin)**2)  # 波束要能够覆盖XZ面内的所有目标
D=5*D   # 可增加系数
Ymin = Ymin-D
Ymax = Ymax+D  # y方向增加 D
T_a=2*D/VR 
TBP_a=BW_a*T_a
if(TBP_a<TBP):
    print('TBP_a is too small') 
    print('Please increase D or decrease VR or decrease Da')
    exit()
#
H1=D/np.tan(Theta_bw) 
if(H1<10*np.linalg.norm(Rmax-Rmin)):
    print('H1 is too small') 
    print('Please increase D')
    exit()
#
ZR=(Zmin+Zmax)/2+H1*np.sin(angle_bw)  # 雷达Z坐标
XR=(Xmin+Xmax)/2-H1*np.cos(angle_bw)  # 雷达X坐标
# HR=D/(sqrt(2)*tan(Theta_bw)) 
# if(HR<10*Zmax)
#     disp('HR is too small') 
#     error('Please increase D')  
# end
#clear D 
###################################################
#RR0=[(Xmin+Xmax)/2-HR Ymin (Zmin+Zmax)/2+HR]  # 雷达初始坐标
RR0=np.array([XR,Ymin,ZR])  # 雷达初始坐标
#
t0=2*np.linalg.norm(np.array([Xmin,Ymin,Zmax])-RR0)/c 
t0=-Tr/2+t0  # 波开始时刻
t1= Tr/2+2*np.linalg.norm(np.array([Xmax,Ymin,Zmin])-RR0)/c  #波结束时刻
if(t1>dt_a-Tr/2):
    print('Echo time > PRI')
    exit()
#
N_r=int(np.ceil((t1-t0)/dt_r))               #距离向像素数
N_a=int(np.ceil((Ymax-Ymin)/(VR*dt_a)))     #方位向像素数
#clear t1
###################################################
with open('mdl.pk','rb') as f:
    NT=load(f) # 目标点个数
    RT=load(f) # 目标点坐标  RT(NT,3)
tr=t0+dt_r*np.arange(0,N_r)
s0=np.zeros([N_r,N_a],dtype=complex)
print('Simulating echo ...')
for ja in range(N_a): #ja=1:N_a
    if(ja%10==0):print(ja,"/",N_a,"...")
    RR=RR0+np.array([0,ja*dt_a*VR,0]) # 当前雷达坐标
    R =RT-np.tile(RR,(NT,1)) # 雷达指向目标点的矢量  R(NT,3)
    w=waveform(R,NT,N_r,d_w,tr,Tr,Kr,F0,Theta_bw)
    s0[:,ja]=s0[:,ja]+np.sum(w,axis=0)
print('Simulate echo end.')
save_s0(s0,dt_r,dt_a,t0,RR0,VR,F0,Kr,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax)
