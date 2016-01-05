#!/usr/bin/python
import sys,os,math
sys.path.append("../")
from ELoss import *

def testH():
    a=ELoss()
    ke=40e3
    E=ke*eV/gE; Mp=131.9*u/gM
    v=math.sqrt(2*E/Mp)
    Zp=1
    n=4.518696e25
    print n,a.calSP(Zp,v,n)*gE/gL/(eV/A)

testH()

