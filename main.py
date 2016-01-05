#!/usr/bin/python
import sys,os,math
from ELoss import *

def testC4C():
    a=ELoss()
    a.readObf("c.orb")
    E=100e3*12*eV/gE; Mp=12*u/gM
    v=math.sqrt(2*E/Mp)
    Zp=4
    dz=0.5*A/gL
    dt=dz/v
    z0=30*A/gL; 

    d0=0.1
    while d0<20.0:
        print d0, a.calELoss(Zp, v, dt, z0, d0*A/gL)*gE/eV
        d0=d0 + 0.1

def testHTi():
    a=ELoss()
    a.readObf("ti.orb")
    ke=1.0
    while ke<=10000:
        E=ke*1e3*eV/gE; Mp=1*u/gM
        v=math.sqrt(2*E/Mp)
        Zp=1
        Mt=47.9*u/gM
        Nt=4.5189/Mt
        print ke*1e3,a.calLDA_SP(Zp,v,Nt)*gE/gL/eV*A
        ke=ke*2

def testCTi():
    a=ELoss()
    a.readObf("ti.orb")
    ke=1.0
    while ke<=10000:
        E=ke*1e3*eV/gE; Mp=12*u/gM
        v=math.sqrt(2*E/Mp)
        Mt=47.9*u/gM
        Nt=4.5189/Mt
        Zp=6
        print ke*1e3,a.calLDA_SP(Zp,v,Nt)*gE/gL/eV*A
        ke=ke*2

