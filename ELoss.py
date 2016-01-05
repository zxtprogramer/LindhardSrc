#!/usr/bin/python
import sys,os,math
from LFun import *

class ELoss:
    def __init__(self):
        self.lfun=LFun()

    def readObf(self, dFile):
        self.dFile=dFile
        self.density=[]
        fp=open(self.dFile,'r')
        for line in fp.readlines():
            items=line.split()
            self.density.append([float(items[0]), float(items[1])])
        fp.close()
        self.dr=self.density[1][0] - self.density[0][0]
        self.length=len(self.density)
        self.rMax=self.density[self.length-1][0]

    def calSP(self, Zp, v, n):
        kf=math.pow(3*PI*PI*n,1.0/3.0)
        vf=hbar_g*kf/me_g
        chi2=e_g*e_g/PI/hbar_g/vf
        chi=math.sqrt(chi2)
        y=math.sqrt(3)/chi*(v/vf)**2
        return self.lfun.calL(y,chi2,1e3)*4*PI*Zp*Zp*e_g**4/me_g/v/v*n

    def calLDA_SP(self, Zp, v, N):
        R0=math.pow(3.0/4.0/PI/N, 1.0/3.0)
        r=0; dr=R0/1.0e2
        sumL=0
        while r<=R0:
            index=int(r/self.dr)
            if index>= self.length:
                break
            else:
                n=float(self.density[index][1])

            kf=math.pow(3*PI*PI*n,1.0/3.0)
            vf=hbar_g*kf/me_g
            chi2=e_g*e_g/PI/hbar_g/vf
            chi=math.sqrt(chi2)
            y=math.sqrt(3)/chi*(v/vf)**2

            sumL=sumL + self.lfun.calL(y,chi2,1e2)*n*4*PI*r*r*dr
            r=r + dr
        return sumL*N*4*PI*Zp*Zp*e_g**4/me_g/v/v
 
        
    def calELoss(self, Zp, v, dt, z0, d0):
        z=-1*abs(z0)
        eLoss=0
        while z<abs(z0):
            dis=math.sqrt(z*z + d0*d0)
            index=int(dis/self.dr)
            if index>= self.length:
                z=z + v*dt
                continue
            else:
                n=float(self.density[index][1])
            eLoss=eLoss + self.calSP(Zp, v, n)*v*dt
            z=z + v*dt
        return eLoss

    def calCS(self, Zp, v, dt, z0):
        r=0; dr=self.rMax/1e2
        cs=0.0
        while r<=self.rMax:
            r=r+dr
            cs=cs + PI*((r+dr)**2 - r**2)*self.calELoss(Zp, v, dt, z0, r)
        return cs
            


