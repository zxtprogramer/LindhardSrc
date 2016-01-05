#!/usr/bin/python
import sys,os,math
sys.path.append("../")
from ELoss import *

def testH():
    a=ELoss()
    ke=100e3
    E=ke*eV/gE; Mp=1*u/gM
    v=math.sqrt(2*E/Mp)
    Zp=1
    n=1e22
    while n<1e32:
        print n,a.calSP(Zp,v,n)*gE/gL/(1000*eV/1e-2)
        n=n*10

testH()

