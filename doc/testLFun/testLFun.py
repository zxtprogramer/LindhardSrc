#!/usr/bin/python

import os,sys,math
sys.path.append("../../")
from LFun import *


a=LFun()

for chi2 in [0.01, 0.02154, 0.046416]:
    print "----",chi2,"----"
    for y in [0.1,0.31623,1.0,1.7783,3.1623,5.6235,10.0,14.678,21.544,31.623,56.235,100.0,316.23,1000.0]:
        L=a.calL(y,chi2,1e3)
        print "y=",y,"L=",L
