#!/usr/bin/python
import sys,os,math
sys.path.append("../")
from ELoss import *

def testH():
    a=ELoss()
    ke=1000e3
    E=ke*eV/gE; Mp=1*u/gM
    v=math.sqrt(2*E/Mp)
    Zp=1
    n=1e22
    while n<1e32:
        print n,a.calSP(Zp,v,n)/n*gE*gL*gL/(eV*1e-24)
        n=n*10

testH()

