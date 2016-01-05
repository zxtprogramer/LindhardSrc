#!/usr/bin/python
import os,sys,math

PI=3.141592654
A=1e-10
cm=1e-2

num=0
density=[]
for iN in range(1,len(sys.argv)):
    fName=sys.argv[iN]
    fp=open(fName, 'r')
    content=fp.readlines()
    if num==0:
        for i in range(0,len(content)):
            line=content[i].split()
            r=float(line[0])
            d=float(line[1])
            density.append([r,d])
    else:
        for i in range(0,len(content)):
            line=content[i].split()
            d=float(line[1])
            density[i][1]=density[i][1] + d
    fp.close()
    num=num +1


density2=[]
sumD=0
dr=density[1][0]-density[0][0]
dr=dr*A/cm
for i in range(0,len(density)-1):
    r=density[i][0]*A/cm
    d=density[i][1]/(A**3)*(cm**3)
    print r,d
    sumD=sumD + 4*PI*r*r*dr*d
print sumD

