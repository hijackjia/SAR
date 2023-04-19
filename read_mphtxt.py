# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 08:25:01 2020

@author: jiash
"""

filename='s.mphtxt'  # 文件名

import numpy as np
from pickle import dump

with open(filename,'r') as f:
    print('reading mphtxt file ...')
    line=f.readline()
    while(not line=='3 # sdim\n'):
        line=f.readline()
    line=f.readline()
    values=line.split()
    NT=int(values[0])
    RT=np.zeros([NT,3],dtype=float)
    for _ in range(3):
        line=f.readline()
    for i in range(NT):
        line=f.readline()
        v=line.split()
        RT[i,:]=np.array([float(v[0]),float(v[1]),float(v[2])])
    print('read mphtxt file end')
with open('mdl.pk','wb') as f:
    dump(NT,f)
    dump(RT,f)
