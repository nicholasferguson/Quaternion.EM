# Coded up equations from article:
# "Physical Space as a Quaternion Structure, I
# Maxwell Equations. A Brief Note" by Jack P.M.

# In this study results of each equation is graphed as a
# 3D prism.

# This suite of equations yield better results when
# all values of quaternions are quasi random to each other....

import math
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import interactive
from plotPrism import plot_prism,plot_prism2  #plotprism.py

def displayResults2(pq,pq2,title):
# Fill in a start prism wire frame, using point and a 3D origin as vertices
# Highlight point as vertex

    p1 = np.array([0.,0.,0.])
    p2 = np.array([pq[1],0.,0.])
    p3 = np.array([0.,pq[2],0.])
    p4 = np.array([0.,0.,pq[3]])

    prism1 = [
        (p1[0],p1[1],p1[2]), (p2[0],p2[1],p2[2]), (p3[0],p3[1],p3[2]), (p4[0],p4[1],p4[2])
    ]

    # Fill in a rotated prism wire frame, using point and a 3D origin as vertices
    # Highlight point as vertex

    p1A = np.array([0.,0.,0.])
    p2A = np.array([pq2[1],0.,0.])
    p3A = np.array([0.,pq2[2],0.])
    p4A = np.array([0.,0.,pq2[3]])

    prism2 = [
        (p1A[0],p1A[1],p1A[2]), (p2A[0],p2A[1],p2A[2]), (p3A[0],p3A[1],p3A[2]), (p4A[0],p4A[1],p4A[2])
    ]

    plot_prism2(prism1,prism2,title)
# Fill in a start prism wire frame, using point and a 3D origin as vertices
# Highlight point as vertex
def displayResults(pq,rq,pq2,title):
        p1 = np.array([0.,0.,0.])
        p2 = np.array([pq[1],0.,0.])
        p3 = np.array([0.,pq[2],0.])
        p4 = np.array([0.,0.,pq[3]])

        prism1 = [
            (p1[0],p1[1],p1[2]), (p2[0],p2[1],p2[2]), (p3[0],p3[1],p3[2]), (p4[0],p4[1],p4[2])
        ]

# Fill in a rotated prism wire frame, using point and a 3D origin as vertices
# Highlight point as vertex

        p1A = np.array([0.,0.,0.])
        p2A = np.array([rq[1],0.,0.])
        p3A = np.array([0.,rq[2],0.])
        p4A = np.array([0.,0.,rq[3]])

        prism2 = [
            (p1A[0],p1A[1],p1A[2]), (p2A[0],p2A[1],p2A[2]), (p3A[0],p3A[1],p3A[2]), (p4A[0],p4A[1],p4A[2])
        ]

        p1B = np.array([0.,0.,0.])
        p2B = np.array([pq2[1],0.,0.])
        p3B = np.array([0.,pq2[2],0.])
        p4B = np.array([0.,0.,pq2[3]])

        prism3 = [
            (p1B[0],p1B[1],p1B[2]), (p2B[0],p2B[1],p2B[2]), (p3B[0],p3B[1],p3B[2]), (p4B[0],p4B[1],p4B[2])
       ]

        plot_prism(prism1,prism2,prism3,title)


def normalize(v, tolerance=0.00001):
    mag2 = sum(n * n for n in v)
    if abs(mag2 - 1.0) > tolerance:
        mag = math.sqrt(mag2)
        v = tuple(n / mag for n in v)
    return np.array(v)

# Equation #3
# a->b  --
def a_Right_b_Eqn03(a,b):
    k= [a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3],
            a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2],
            a[0]*b[2]-a[1]*b[3]+a[2]*b[0]+a[3]*b[1],
            a[0]*b[3]+a[1]*b[2]-a[2]*b[1]+a[3]*b[0]]
    return k

# Equation #4
# b <- a  --
def b_Left_a_Eqn04(b,a):
    k= [a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3],
            a[0]*b[1]+a[1]*b[0]-a[2]*b[3]+a[3]*b[2],
            a[0]*b[2]+a[1]*b[3]+a[2]*b[0]-a[3]*b[1],
            a[0]*b[3]-a[1]*b[2]+a[2]*b[1]+a[3]*b[0]]
    return k

# Equation #5
# {a, b} = (1/2)(a → b + b ← a) (5)

def symmetric_Eqn05x(a,b):
    R = a_Right_b_Eqn03(a,b)
    L = b_Left_a_Eqn04(b,a)
    k= [.5*(R[0]+L[0]),
            .5*(R[1]+L[1]),
            .5*(R[2]+L[2]),
            .5*(R[3]+L[3])]
    return k

# Equation #5
# version as per Jacks article
def symmetric_Eqn05(a,b):
    k= [.5*(a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3]),
            .5*(a[0]*b[1]+a[1]*b[0]),
            .5*(a[0]*b[2]+a[2]*b[0]),
            .5*(a[0]*b[0]+a[3]*b[0])]
    return k

# Equation #6
# [a, b] = (1/2)(a → b − b ← a)
def antisymmetric_Eqn06x(a,b):
    A = a
    dr = b
    print("antisymmetric_Eqn06(a,b) a ",a)
    print("antisymmetric_Eqn06(a,b) b ",b)
    R = a_Right_b_Eqn03(A,dr)
    L = b_Left_a_Eqn04(dr,A)
    Z = [.5*(R[0]-L[0]),
            .5*(R[1]-L[1]),
            .5*(R[2]-L[2]),
            .5*(R[3]-L[3])]
    return Z

# Equation #6
# version as per Jacks article
# [a, b] = (1/2)(a → b − b ← a) (6)
def antisymmetric_Eqn06(a,b):

    print("antisymmetric_Eqn06(a,b) a ",a)
    print("antisymmetric_Eqn06(a,b) b ",b)
    Z = [0,
            .5*(a[2]*b[3]-a[3]*b[2]),
            .5*(a[3]*b[1]-a[1]*b[3]),
            .5*(a[1]*b[2]-a[2]*b[1])]
    return Z

# Equation #9
# electric field is the negative symmetric
# derivative of the potential

# E = −{d/dr, A} = −(1/2)(d/dr → A + A ← d/dr) (9)
def fE_Eqn09(a,b):
    # r = normalize(rq)
     A=a
     dr=b
     R = a_Right_b_Eqn03(dr,A)
     L = b_Left_a_Eqn04(A,dr)
     return [-1*.5*(R[0]+L[0]),
             -1*.5*(R[1]+L[1]),
             -1*.5*(R[2]+L[2]),
             -1*.5*(R[3]+L[3])]

# Equation #10
# magnetic field is the positive antisymmetric_Eqn06
# derivative of the potential

# B = +[d/dr, A] = +(1/2)(d/dr → A − A ← d/dr) (10)
def fB_Eqn10(a,b):
    #r = normalize(rq)
     A=a
     dr = b
     R = a_Right_b_Eqn03(dr,A)
     L = b_Left_a_Eqn04(A,dr)
     return [.5*(R[0]-L[0]),
             .5*(R[1]-L[1]),
             .5*(R[2]-L[2]),
             .5*(R[3]-L[3])]

# Equation #13
#  [d/dr, B] = +{d/dr, E} (13)
def Equation13(a,b):
    d_dr = b
    A = a
    BB=antisymmetric_Eqn06(d_dr,fB_Eqn10(A,d_dr))
    EE=symmetric_Eqn05(d_dr,fE_Eqn09(A,d_dr))
    print("Eqn 13: antisymmetric_Eqn06 [d/dr,B]",BB)
    print("Eqn 13: symmetric_Eqn05 +{d/dr,E}",EE)
    ret = [0,0,0,0]
    EE[0] = abs(EE[0])
    EE[1] = abs(EE[1])
    EE[2] = abs(EE[2])
    EE[3] = abs(EE[3])
    ret[0] = BB[0] - EE[0]
    ret[1] = BB[1] - EE[1]
    ret[2] = BB[2] - EE[2]
    ret[3] = BB[3] - EE[3]
    #r = normalize(rq)
    #a_conj = [a[0],-1*a[1],-1*a[2],-1*a[3]]
    #b_conj = [b[0],-1*b[1],-1*b[2],-1*b[3]]
    return BB,EE,ret

# Equation #14
# [d/dr, E] = −{d/dr, B} (14)
def Equation14(a,b):
    d_dr = b
    A = a
    BB=antisymmetric_Eqn06(d_dr,fE_Eqn09(A,d_dr))
    EE=symmetric_Eqn05(d_dr,fB_Eqn10(A,d_dr))
    print("Eqn 13: antisymmetric_Eqn06 [d/dr,B]",BB)
    print("Eqn 13: symmetric_Eqn05 +{d/dr,E}",EE)
    ret = [0,0,0,0]
    ret[0] = BB[0] -  (-1*EE[0])
    ret[1] = BB[1] -  (-1*EE[1])
    ret[2] = BB[2] -  (-1*EE[2])
    ret[3] = BB[3] -  (-1*EE[3])
    #r = normalize(rq)
    #a_conj = [a[0],-1*a[1],-1*a[2],-1*a[3]]
    #b_conj = [b[0],-1*b[1],-1*b[2],-1*b[3]]
    return BB,EE,ret

# Equation #22
# d/dr → (d/dr → A) + (A ← d/dr) ← d/dr = 0 (22)
def Equation22(a,b):
    A = a
    dr = b
    e23=a_Right_b_Eqn03(dr,a_Right_b_Eqn03(dr,A))
    e24=b_Left_a_Eqn04(b_Left_a_Eqn04(A,dr),dr)
    print("eq22 e23",e23)
    print("eq22 e24",e24)
    e25=[0,0,0,0]
    e25[0]=e23[0]+e24[0]
    e25[1]=e23[1]+e24[1]
    e25[2]=e23[2]+e24[2]
    e25[3]=e23[3]+e24[3]

    return e25;

# Equation #23
# d/dr → (A ← d/dr) − (d/dr → A) ← d/dr = 0 (23)
def Equation23(a,b):
   # A = normalize(a)
   # dr = normalize(b)
    A = a
    dr = b
    e23=a_Right_b_Eqn03(dr,b_Left_a_Eqn04(A,dr))
    print("eqn 23 e23",e23)
    e24=b_Left_a_Eqn04(a_Right_b_Eqn03(dr,A),dr)
    print("eqn 23 e24",e24)
    e25=[0,0,0,0]
    e25[0]=e23[0]-e24[0]
    e25[1]=e23[1]-e24[1]
    e25[2]=e23[2]-e24[2]
    e25[3]=e23[3]-e24[3]

    return e25;

# Equation #24
#  d/dr → (d/dr → A) + (A ← d/dr) ← d/dr = 8πJ (24)
def Equation24(a,b):
    A = a
    dr = b
    e23=a_Right_b_Eqn03(dr,a_Right_b_Eqn03(dr,A))
    e24=b_Left_a_Eqn04(b_Left_a_Eqn04(A,dr),dr)
    print("eq24 e23",e23)
    print("eq24 e24",e24)
    e25=[0,0,0,0]
    e25[0]=e23[0]+e24[0]
    e25[1]=e23[1]+e24[1]
    e25[2]=e23[2]+e24[2]
    e25[3]=e23[3]+e24[3]

    return e25;

def main():

    A = [.01, .03, .06, .14]  # Electric potential.
    d_dr = [.001,.0025,.0011,.07]  # d/dr.
    interactive(True)

    #r = normalize(rq)
    AB = a_Right_b_Eqn03(A,d_dr)
    print("a_Right_b_Eqn03 AB",AB)
    displayResults(A,d_dr,AB,"a_Right_b_Eqn03")


    BA = b_Left_a_Eqn04(A,d_dr)
    print("b_Left_a_Eqn04 ",BA)
    displayResults(A,d_dr,BA,"b_Left_a_Eqn04")

    aAB = antisymmetric_Eqn06(A,d_dr)
    print("antisymmetric_Eqn06 ",aAB)
    displayResults(A,d_dr,aAB,"antisymmetric_Eqn06")

    sAB = symmetric_Eqn05(A,d_dr)
    print("symmetric_Eqn05 ",sAB)
    displayResults(A,d_dr,sAB,"symmetric_Eqn05")

    B = fB_Eqn10(A,d_dr)
    print("B :",B)
    displayResults(A,d_dr,B,"B Eqn: 10")

    E = fE_Eqn09(A,d_dr)
    print("E :",E)
    displayResults(A,d_dr,E,"E Eqn: 09")

    BB,EE,ret = Equation13(A,d_dr)

    interactive(False)

    displayResults2(BB,EE,"[d/dr,B] +{d/dr,E}")



    print("Equation22 \n",Equation22(A,d_dr))

    print("Equation23 \n",Equation23(A,d_dr))

    print("Equation24 \n",Equation24(A,d_dr))



if __name__ == "__main__":
    main()

