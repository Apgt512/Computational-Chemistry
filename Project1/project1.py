# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 23:26:21 2016

@author: apgt
"""
import numpy as np

file=open("coords.dat","r")
output = open('output.txt','w')
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
def dist(i,j):
        R=np.sqrt((coords[i][1]-coords[j][1])**2+(coords[i][2]-\
        coords[j][2])**2+(coords[i][3]-coords[j][3])**2)
        return R
        
#Calculate angles between atoms i,j,k. The central atom is j
def angle(i,j,k):
    e1=[]
    e2=[]
    for c in range(1,4):
        e1.append(-(coords[j][c]-coords[i][c])/dist(i,j))
        e2.append(-(coords[j][c]-coords[k][c])/dist(j,k))
    phi=np.arccos(np.dot(e1,e2))
    return phi

#Calculate angle between atom i and plane j-k-l
def oop_ang(i,j,k,l):
    e1,e2,e3=[],[],[]
    for c in range(1,4):
        e1.append(-(coords[k][c]-coords[j][c])/dist(k,j))
        e2.append(-(coords[k][c]-coords[l][c])/dist(k,l))
        e3.append(-(coords[k][c]-coords[i][c])/dist(k,i))
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
        e1.append(-(coords[i][c]-coords[j][c])/dist(i,j))
        e2.append(-(coords[j][c]-coords[k][c])/dist(j,k))
        e3.append(-(coords[k][c]-coords[l][c])/dist(k,l))
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
    for i in range(len(coords)):
        for j in range(3):
            coords[i][j+1] = coords[i][j+1] - center_of_mass[j]
    return center_of_mass

#Calculate the moment of inertia
def inertia(xyz):
    I_xx, I_yy, I_zz = 0, 0, 0 #Diagonal
    I_xy, I_xz, I_yz = 0, 0, 0 #Off diagonal
    for i in range(len(xyz)):
        mass = atomic_mass[xyz[i][0]]
        I_xx += mass * (xyz[i][2]**2 + xyz[i][3]**2)
        I_yy += mass * (xyz[i][1]**2 + xyz[i][3]**2)
        I_zz += mass * (xyz[i][1]**2 + xyz[i][2]**2)
        I_xy += mass * xyz[i][1] * xyz[i][2]
        I_xz += mass * xyz[i][1] * xyz[i][3]
        I_yz += mass * xyz[i][2] * xyz[i][3]
    I = np.array([I_xx,I_xy,I_xz,I_xy,I_yy,I_yz,I_xz,I_yz,I_zz]).reshape(3,3)
    output.write('\nMoment of inertia tensor:\n{}\n'.format(I))
    w, v = np.linalg.eig(I)
    w = np.sort(w)
    output.write('\nPrincipal moments of inertia:\n{}\n'.format(w))
    if np.abs(w[0]-w[1])<1E-4 and np.abs(w[1]-w[2])<1E-4:
        output.write('\nMolecule is a spherical top\n')
    elif np.abs(w[0]-w[1])<1E-4 and np.abs(w[1]-w[2])>1E-4:
        output.write('\nMolecule is an oblate symmetric top\n')
    elif np.abs(w[0]-w[1])>1E-4 and np.abs(w[1]-w[2])<1E-4:
        output.write('\nMolecule is a prolate symmetric top\n')
    else:
        output.write('\nMolecule is an asymmetric top\n')
    pi = np.pi
    conv = 6.6260755E-34/(8*pi**2)
    conv /= 1.6605402E-27 * 0.529177249E-10 * 0.529177249E-10
    conv *= 1E-6
    #Rotational constants
    A, B, C = conv/w[0], conv/w[1], conv/w[2]
    output.write('\nRotational constants (MHz):\nA = {:.2f}\tB = {:.2f}\tC = {:.2f}'.format(A,B,C))

##########################Print the results############################
#Print distance
output.write('Interatomic distances:\n')
for i in range(len(coords)):
    for j in range(i):
        output.write('{} {} {:.6f}\n'.format(i,j,dist(i,j)))
print('1/6')

#Print bond angles
output.write('\nBond angles:\n')
for i in range(len(coords)):
    for j in range(i):
        for k in range(j):
            if dist(i,j)<4.0 and dist(j,k)<4.0:
                output.write('{}-{}-{} {:.4f}\n'\
                .format(i,j,k,angle(i,j,k)*(180.0/np.arccos(-1.0))))
print('2/6')

#Print out of plane angles
output.write('\nOut of plane angles:\n')
for i in range(len(coords)):
    for k in range(len(coords)):
        for j in range(len(coords)):
            for l in range(j):
                if (i!=j and i!=k and i!=l and j!=k and k!=l and\
                dist(i,k)<4.0 and dist(k,j)<4.0 and dist(k,l)<4.0):
                    output.write('{}-{}-{}-{} {:.4f}\n'\
                    .format(i,j,k,l,oop_ang(i,j,k,l)*(180.0/np.arccos(-1.0))))
print('3/6')

#Print torsional angles
output.write('\nTorsional angles:\n')
for i in range(len(coords)):
    for j in range(i):
        for k in range(j):
            for l in range(k):
                if (dist(i,j)<4.0 and dist(j,k)<4.0 and dist(k,l)<4.0):
                    output.write('{}-{}-{}-{} {:.4f}\n'\
                    .format(i,j,k,l,dihedral(i,j,k,l)*(180.0/np.arccos(-1.0))))
print('4/6')

#Print center of mass
output.write('\nCenter of mass:\n')
output.write('{}\n'.format(com(coords)))
print('5/6')

#Print moment of inertia
inertia(coords)
print('6/6')

output.close()
file.close()