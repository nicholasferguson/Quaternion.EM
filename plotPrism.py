import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

def new_subplot(self):
                nsubplots = len(self.__subplots) + 1

                for i, subplot in enumerate(self.__subplots):
                        subplot.change_geometry(nsubplots, 1, i + 1)

                subplot = self.figure.add_subplot(nsubplots, 1, nsubplots)
                subplot.grid(True)

                self.__subplots.append(subplot)
                self.__subplot = subplot

def plot_prism(prism_definition1, prism_definition2,prism_definition3,title):
    prism_definition_array1 = [
        np.array(list(item))
        for item in prism_definition1
    ]
    prism_definition_array2 = [
        np.array(list(item))
        for item in prism_definition2
    ]
    prism_definition_array3 = [
        np.array(list(item))
        for item in prism_definition3
    ]
    pt1x = (prism_definition_array1[1][0])
    pt1y = (prism_definition_array1[2][1])
    pt1z = (prism_definition_array1[3][2])

    pt2x = (prism_definition_array2[1][0])
    pt2y = (prism_definition_array2[2][1])
    pt2z = (prism_definition_array2[3][2])

    pt3x = (prism_definition_array3[1][0])
    pt3y = (prism_definition_array3[2][1])
    pt3z = (prism_definition_array3[3][2])
#================================
    points1 = []
    points1 += prism_definition_array1
    vectors1 = [
        prism_definition_array1[1] - prism_definition_array1[0],
        prism_definition_array1[2] - prism_definition_array1[0],
        prism_definition_array1[3] - prism_definition_array1[0]
    ]

    points1 += [prism_definition_array1[0] + vectors1[0] + vectors1[1]]
    points1 += [prism_definition_array1[0] + vectors1[0] + vectors1[2]]
    points1 += [prism_definition_array1[0] + vectors1[1] + vectors1[2]]
    points1 += [prism_definition_array1[0] + vectors1[0] + vectors1[1] + vectors1[2]]

    points1 = np.array(points1)

    edges1 = [
        [points1[0], points1[3], points1[5], points1[1]],
        [points1[1], points1[5], points1[7], points1[4]],
        [points1[4], points1[2], points1[6], points1[7]],
        [points1[2], points1[6], points1[3], points1[0]],
        [points1[0], points1[2], points1[4], points1[1]],
        [points1[3], points1[6], points1[7], points1[5]]
    ]

#================================
    points2 = []
    points2 += prism_definition_array2
    vectors2 = [
        prism_definition_array2[1] - prism_definition_array2[0],
        prism_definition_array2[2] - prism_definition_array2[0],
        prism_definition_array2[3] - prism_definition_array2[0]
    ]

    points2 += [prism_definition_array2[0] + vectors2[0] + vectors2[1]]
    points2 += [prism_definition_array2[0] + vectors2[0] + vectors2[2]]
    points2 += [prism_definition_array2[0] + vectors2[1] + vectors2[2]]
    points2 += [prism_definition_array2[0] + vectors2[0] + vectors2[1] + vectors2[2]]

    points2 = np.array(points2)

    edges2 = [
        [points2[0], points2[3], points2[5], points2[1]],
        [points2[1], points2[5], points2[7], points2[4]],
        [points2[4], points2[2], points2[6], points2[7]],
        [points2[2], points2[6], points2[3], points2[0]],
        [points2[0], points2[2], points2[4], points2[1]],
        [points2[3], points2[6], points2[7], points2[5]]
    ]

#================================
    points3 = []
    points3 += prism_definition_array3
    vectors3 = [
        prism_definition_array3[1] - prism_definition_array3[0],
        prism_definition_array3[2] - prism_definition_array3[0],
        prism_definition_array3[3] - prism_definition_array3[0]
    ]

    points3 += [prism_definition_array3[0] + vectors3[0] + vectors3[1]]
    points3 += [prism_definition_array3[0] + vectors3[0] + vectors3[2]]
    points3 += [prism_definition_array3[0] + vectors3[1] + vectors3[2]]
    points3 += [prism_definition_array3[0] + vectors3[0] + vectors3[1] + vectors3[2]]

    points3 = np.array(points3)
    edges3 = [
        [points3[0], points3[3], points3[5], points3[1]],
        [points3[1], points3[5], points3[7], points3[4]],
        [points3[4], points3[2], points3[6], points3[7]],
        [points3[2], points3[6], points3[3], points3[0]],
        [points3[0], points3[2], points3[4], points3[1]],
        [points3[3], points3[6], points3[7], points3[5]]
    ]
#===============================
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    faces1 = Poly3DCollection(edges1, linewidths=1, edgecolors='k')
    faces1.set_facecolor((0,1,1,0.1))  #greem dot
    ax.add_collection3d(faces1)

    faces2 = Poly3DCollection(edges2, linewidths=1, edgecolors='k')
    faces2.set_facecolor((0,0,1,0.1)) #blue  dot
    ax.add_collection3d(faces2)

    faces3 = Poly3DCollection(edges3, linewidths=1, edgecolors='k')
    faces3.set_facecolor(((0.8, 0.6, 0.5,0.1))) # red  dot
    ax.add_collection3d(faces3)

    labelW='test'
    if title == 'a_Right_b_Eqn03':
        labelW='a->b'
        labelA ='a=A'
        labelB ='b=d_dr'
    elif title == 'b_Left_a_Eqn04':
        labelW='b<-a'
        labelA ='a=A'
        labelB ='b=d_dr'
    elif title == 'symmetric_Eqn05':
        labelW = '1/2[a->b + b<-a] {exp}'
        labelA ='a=A'
        labelB ='b=d_dr'
    elif title == 'antisymmetric_Eqn06':
        labelW = '1/2[a->b - b<-a] [exp]'
        labelA ='a=A'
        labelB ='b=d_dr'
    elif title == 'E Eqn: 09':
        labelW = '  -1* 1/2[a->b + b<-a] {exp}'
        labelA ='a=A'
        labelB ='b=d_dr'
    elif title == 'B Eqn: 10':
        labelW = '1/2[a->b - b<-a] [exp]'
        labelA ='a=A'
        labelB ='b=d_dr'
    elif title =='[d/dr,B] +{d/dr,E}':
        title = '+{d/dr,E}  [d/dr,B]'
        labelW = '[d/dr,B] +{d/dr,E})'
        labelA ='[d/dr,B] antisymmetric(d/dr,B=[d/dr,A])'
        labelB ='+{d/dr,E} symmetric(d/dr,E=-{d/dr,A})'




    if title == '[d/dr,B]-+{d/dr,E}':
        ax.plot(pt1x, pt1y, pt1z, markerfacecolor='g', markeredgecolor='k', marker='o', markersize=10, alpha=0.6,label=labelA)
        ax.plot(pt2x, pt2y, pt2z, markerfacecolor='b', markeredgecolor='k', marker='o', markersize=10, alpha=0.6,label=labelB)
    else:
        ax.plot(pt1x, pt1y, pt1z, markerfacecolor='g', markeredgecolor='k', marker='o', markersize=10, alpha=0.6,label=labelA)
        ax.plot(pt2x, pt2y, pt2z, markerfacecolor='b', markeredgecolor='k', marker='o', markersize=10, alpha=0.6,label=labelB)
        ax.plot(pt3x, pt3y, pt3z, markerfacecolor='r', markeredgecolor='k', marker='o', markersize=10, alpha=0.6,label=labelW)

   # maxX = max(pt1x,pt2x,pt3x)
   # maxY= max(pt1y,pt2y,pt3y)
   # maxY= max(pt1z,pt2z,pt3z)
    maxA = max(abs(pt1x),abs(pt2x),abs(pt3x),abs(pt1y),abs(pt2y),abs(pt3y),abs(pt1z),abs(pt2z),abs(pt3z))

    # Plot the points themselves to force the scaling of the axes
  #  ax.scatter(points1[:,0], points1[:,1], points1[:,2], s=1)

 #   ax.set_aspect('auto')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Setting the axes properties
    maxA = maxA*2
    ax.set_xlim3d([-maxA, maxA])
    ax.set_ylim3d([-maxA, maxA])
    ax.set_zlim3d([-maxA, maxA])
  #  return ax
    plt.title(title)
    plt.legend(loc='upper right')

    plt.show()


def plot_prism2(prism_definition1, prism_definition2,title):
    prism_definition_array1 = [
        np.array(list(item))
        for item in prism_definition1
    ]
    prism_definition_array2 = [
        np.array(list(item))
        for item in prism_definition2
    ]
    pt1x = (prism_definition_array1[1][0])
    pt1y = (prism_definition_array1[2][1])
    pt1z = (prism_definition_array1[3][2])

    pt2x = (prism_definition_array2[1][0])
    pt2y = (prism_definition_array2[2][1])
    pt2z = (prism_definition_array2[3][2])

#================================
    points1 = []
    points1 += prism_definition_array1
    vectors1 = [
        prism_definition_array1[1] - prism_definition_array1[0],
        prism_definition_array1[2] - prism_definition_array1[0],
        prism_definition_array1[3] - prism_definition_array1[0]
    ]

    points1 += [prism_definition_array1[0] + vectors1[0] + vectors1[1]]
    points1 += [prism_definition_array1[0] + vectors1[0] + vectors1[2]]
    points1 += [prism_definition_array1[0] + vectors1[1] + vectors1[2]]
    points1 += [prism_definition_array1[0] + vectors1[0] + vectors1[1] + vectors1[2]]

    points1 = np.array(points1)

    edges1 = [
        [points1[0], points1[3], points1[5], points1[1]],
        [points1[1], points1[5], points1[7], points1[4]],
        [points1[4], points1[2], points1[6], points1[7]],
        [points1[2], points1[6], points1[3], points1[0]],
        [points1[0], points1[2], points1[4], points1[1]],
        [points1[3], points1[6], points1[7], points1[5]]
    ]

#================================
    points2 = []
    points2 += prism_definition_array2
    vectors2 = [
        prism_definition_array2[1] - prism_definition_array2[0],
        prism_definition_array2[2] - prism_definition_array2[0],
        prism_definition_array2[3] - prism_definition_array2[0]
    ]

    points2 += [prism_definition_array2[0] + vectors2[0] + vectors2[1]]
    points2 += [prism_definition_array2[0] + vectors2[0] + vectors2[2]]
    points2 += [prism_definition_array2[0] + vectors2[1] + vectors2[2]]
    points2 += [prism_definition_array2[0] + vectors2[0] + vectors2[1] + vectors2[2]]

    points2 = np.array(points2)

    edges2 = [
        [points2[0], points2[3], points2[5], points2[1]],
        [points2[1], points2[5], points2[7], points2[4]],
        [points2[4], points2[2], points2[6], points2[7]],
        [points2[2], points2[6], points2[3], points2[0]],
        [points2[0], points2[2], points2[4], points2[1]],
        [points2[3], points2[6], points2[7], points2[5]]
    ]


#===============================
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    faces1 = Poly3DCollection(edges1, linewidths=1, edgecolors='k')
    faces1.set_facecolor((0,1,1,0.1))  #greem dot
    ax.add_collection3d(faces1)

    faces2 = Poly3DCollection(edges2, linewidths=1, edgecolors='k')
    faces2.set_facecolor((0,0,1,0.1)) #blue  dot
    ax.add_collection3d(faces2)


    labelW='test'
    if title == 'a_Right_b_Eqn03':
        labelW='a->b'
        labelA ='a=A'
        labelB ='b=d_dr'
    elif title == 'b_Left_a_Eqn04':
        labelW='b<-a'
        labelA ='a=A'
        labelB ='b=d_dr'
    elif title == 'symmetric_Eqn05':
        labelW = '1/2[a->b + b<-a] {exp}'
        labelA ='a=A'
        labelB ='b=d_dr'
    elif title == 'antisymmetric_Eqn06':
        labelW = '1/2[a->b - b<-a] [exp]'
        labelA ='a=A'
        labelB ='b=d_dr'
    elif title == 'E Eqn: 09':
        labelW = '  -1* 1/2[a->b + b<-a] {exp}'
        labelA ='a=A'
        labelB ='b=d_dr'
    elif title == 'B Eqn: 10':
        labelW = '1/2[a->b - b<-a] [exp]'
        labelA ='a=A'
        labelB ='b=d_dr'
    elif title =='[d/dr,B] +{d/dr,E}':
        title = '+{d/dr,E}  [d/dr,B]'
        labelW = '[d/dr,B] +{d/dr,E})'
        labelA ='[d/dr,B] antisymmetric(d/dr,B=[d/dr,A])'
        labelB ='+{d/dr,E} symmetric(d/dr,E=-{d/dr,A})'





        ax.plot(pt1x, pt1y, pt1z, markerfacecolor='g', markeredgecolor='k', marker='o', markersize=10, alpha=0.6,label=labelA)
        ax.plot(pt2x, pt2y, pt2z, markerfacecolor='b', markeredgecolor='k', marker='o', markersize=10, alpha=0.6,label=labelB)
    # maxX = max(pt1x,pt2x,pt3x)
   # maxY= max(pt1y,pt2y,pt3y)
   # maxY= max(pt1z,pt2z,pt3z)
    maxA = max(abs(pt1x),abs(pt2x),abs(pt1y),abs(pt2y),abs(pt1z),abs(pt2z))

    # Plot the points themselves to force the scaling of the axes
  #  ax.scatter(points1[:,0], points1[:,1], points1[:,2], s=1)

 #   ax.set_aspect('auto')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Setting the axes properties
    maxA = maxA*2
    ax.set_xlim3d([-maxA, maxA])
    ax.set_ylim3d([-maxA, maxA])
    ax.set_zlim3d([-maxA, maxA])
  #  return ax
    plt.title(title)
    plt.legend(loc='upper right')

    plt.show()