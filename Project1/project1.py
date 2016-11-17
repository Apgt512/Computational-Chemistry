# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 23:26:21 2016

@author: apgt
"""
import numpy as np

file=open("coords.dat","r")
read_coords=[]
coords=[]
num_atoms=0

#Atomic mass
atomic_mass={1:1.008,6:12.010,8:15.999,79:196.966}

#To print output of angles
#*(180.0/np.arccos(-1.0))

#Read file
for line in file:
    try:
        num_atoms=int(line)
    except:
        read_coords.append(line.split())


#Convert coordinates from file to float
for i in range(num_atoms):
    lst=[]
    lst=[float(j) for j in read_coords[i]]
    coords.append(lst)


#Calculate distance between atom i and j
def distance(i,j):
        R=np.sqrt((coords[i][1]-coords[j][1])**2+(coords[i][2]-\
        coords[j][2])**2+(coords[i][3]-coords[j][3])**2)
        return R
        
#Calculate angles between atoms i,j,k. The central atom is j
def angle(i,j,k):
    e1=[]
    e2=[]
    for c in range(1,4):
        e1.append(-(coords[j][c]-coords[i][c])/distance(i,j))
        e2.append(-(coords[j][c]-coords[k][c])/distance(j,k))
    phi=np.arccos(np.dot(e1,e2))
    return phi

#Calculate angle between atom i and plane j-k-l
def oop_ang(i,j,k,l):
    e1,e2,e3=[],[],[]
    for c in range(1,4):
        e1.append(-(coords[k][c]-coords[j][c])/distance(k,j))
        e2.append(-(coords[k][c]-coords[l][c])/distance(k,l))
        e3.append(-(coords[k][c]-coords[i][c])/distance(k,i))
    theta=np.cross(e1,e2)
    theta=theta/(np.sin(angle(j,k,l)))
    theta=np.dot(theta,e3)
    if theta<-1.0:
        theta=np.arcsin(-1.0)
    elif theta>1.0:
        theta=np.arcsin(1.0)
    else:
        theta=np.arcsin(theta)
    return theta

#Calculate dihedral angles
def dihedral(i,j,k,l):
    e1,e2,e3 = [],[],[]
    for c in range(1,4):
        e1.append(-(coords[i][c]-coords[j][c])/distance(i,j))
        e2.append(-(coords[j][c]-coords[k][c])/distance(j,k))
        e3.append(-(coords[k][c]-coords[l][c])/distance(k,l))
    ang1=np.sin(angle(i,j,k))
    ang2=np.sin(angle(j,k,l))
    tau=np.dot(np.cross(e1,e2),np.cross(e2,e3))
    tau=tau/(ang1*ang2)
    if tau<-1.0:
        tau=np.arccos(-1.0)
    elif tau>1.0:
        tau=np.arccos(1.0)
    else:
        tau=np.arccos(tau)
    return tau

#Calculate center of mass
def com(xyz):
    total_mass=0
    for i in range(len(xyz)):
        total_mass+=atomic_mass[xyz[i][0]]
    center_of_mass=[0,0,0]
    for i in range(1,4):
        for j in range(len(xyz)):
            center_of_mass[i-1]+=xyz[j][i]*atomic_mass[xyz[j][0]]/total_mass
    return center_of_mass