#!/usr/bin/python

import os,sys,math
from scipy import integrate

A=1e-10
aB=0.529177*A
PI=3.14159265354
h=6.62606957e-34
hbar=h/2.0/PI
me=9.10938291e-31
e=1.60217656e-19
eps0=8.854187e-12
u=1.6605655e-27
eV=e


gE=1e-7
gL=1e-2
gF=1e-5
gM=1e-3
gC=3.33564096e-10
gV=gL/1

aB_g=aB/gL;
h_g=h/gE
hbar_g=h_g/2/PI
me_g=me/gM
e_g=e/gC


class epsilon:
    'use gauss unit'
    def __init__(self):
        pass

    def alpha1(self,z,u):
        a=z-u; b=z+u
        p1=0.5
        if abs(a)==1.0:
            p2=0
        else:
            p2=1.0/8.0/z*(1.0-a*a)*math.log(abs((a+1.0)/(a-1.0))) 
            
        if abs(b)==1.0:
            p3=0
        else:
            p3=1.0/8.0/z*(1-b*b)*math.log(abs((b+1.0)/(b-1.0)))

        return p1+p2+p3

    def alpha2(self,z,u):
        a=z-u; b=z+u
        if b<=1:
            return PI/2.0*u
        if b>1 and abs(a)<1:
            return PI/8.0/z*(1-a*a)
        if abs(a)>=1:
            return 0

    def eReal(self,k,omega,kf):
        vf=hbar_g*kf/me_g
        z=k/2.0/kf; u=omega/k/vf;
        chi2=e_g*e_g/PI/hbar_g/vf
        return 1+chi2/z/z*self.alpha1(z,u)
        
        
    def eImg(self,k,omega,kf):
        vf=hbar_g*kf/me_g
        z=k/2.0/kf; u=omega/k/vf;
        chi2=e_g*e_g/PI/hbar_g/vf
        return chi2/z/z*self.alpha2(z,u)

    def fun1(self,k,omega,kf):
        r=self.eReal(k,omega,kf)
        i=self.eImg(k,omega,kf)
        return -1*i/(r*r + i*i)


    def testA(self):
        z=0; u=0
        vf=e_g*e_g/0.1/PI/hbar_g; kf=vf*me_g/hbar_g
        #kf=10/PI/aB_g; vf=hbar_g*kf/me_g

        while u<3:
            u=u+0.01
            z=0
            while z<3:
                z=z+0.01
                k=z*2.0*kf; omega=k*u*vf;
                print z,u,-1*z*self.fun1(k,omega,kf)*1e3

class LFun:
    def __init__(self):
        self.eps=epsilon()

    def ucFun(self,chi2, u):
        return (u-1)**2 + chi2/2.0*(1 - u*math.log(u/(u-1)))

    def calZeroPoint(self,fun,u0,u1,pre):
        f0=fun(u0); f1=fun(u1)
        du=u1-u0
        while du>pre and f0*f1<0:
            f0=fun(u0); f1=fun(u1)
            u=(u0+u1)/2.0
            f=fun(u)
            if f*f0<=0:
                u1=u
            if f*f1<=0:
                u0=u
            du=u1-u0
        return (u0+u1)/2.0

    def calL(self,y,chi2,divNum):
        vf=e_g*e_g/chi2/PI/hbar_g
        chi=math.sqrt(chi2)
        v=math.sqrt(y*chi/math.sqrt(3))*vf
        kf=vf*me_g/hbar_g

        zMax=v/vf+1

#        uc=math.sqrt(0.25+chi/math.sqrt(3))+0.5
#        if v/vf>uc:
#            return math.log(y)-math.pow(3,1.5)/5/chi/y - 9.0/14.0/chi/chi/y/y


        f1=lambda u: self.ucFun(chi2,u)
        pre=1e-3
        uc=self.calZeroPoint(f1,1.0+pre**2,10.0,pre)
        zc=uc-1.0

        if v/vf>uc:
            return math.log(y)-math.pow(3,1.5)/5/chi/y - 9.0/14.0/chi/chi/y/y


        du=v/vf/divNum
        dz=zMax/divNum
        u=0.0; z=0.0;
        sumZ=0.0; sumU=0.0
        while u<=(v/vf):
            u=u + du
            sumZ=0
            z=0
            #print u/(v/vf)*100,"%"
            while z<=zMax:
                z=z + dz
                if abs(u-z)>=1:
                    continue
                k=2*z*kf; omega=u*k*vf
                sumZ=sumZ + z*self.eps.fun1(k,omega,kf)*dz
            sumU=sumU + u*du*sumZ

        return -6*sumU/PI/chi2

