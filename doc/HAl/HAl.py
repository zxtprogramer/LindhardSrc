#!/usr/bin/python
import sys,os,math 
sys.path.append("../../")
from ELoss import *

def testHAl():
    a=ELoss()
    a.readObf("al.orb_fit")
    ke=10000e3
    E=ke*eV/gE; Mp=1.0*u/gM
    v=math.sqrt(2*E/Mp)
    Zp=1.0
    i=0
    while i<len(a.density):
        d=float(a.density[i][0])
        n=float(a.density[i][1])
        d0=d*gL/A
        nm=A*10.0
        sp=a.calSP(Zp,v,n)
        L=a.calLFun(Zp,v,n)
        #print d0,n*4*PI*d**2,L,sp,sp*4*PI*d*d*gE*gL/eV/A
        print d0,L,L*n*4*PI*d*d/gL*A
        i=i+1
    
testHAl()

