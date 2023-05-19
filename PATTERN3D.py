##3D PATTERN

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d

##NODES AND CONNECTIONS LATTICE DEFINITION
nodes = np.array([[0,0,0],[0,1,0],[1,0,0],[1,1,0],[0.5,1,0],[0.5,0,0],[1.5,0.5,0],[-0.5,0.5,0],[0.5,0.5,0],[3.333333333333333148e-01,8.333333333333333703e-01,1],[6.666666666666666297e-01,8.333333333333333703e-01,1],[3.333333333333333148e-01,1.666666666666666574e-01,1],[6.666666666666666297e-01,1.666666666666666574e-01,1],[0.000000000000000000e+00,6.666666666666666297e-01,1],[0.000000000000000000e+00,3.333333333333333148e-01,1],[1.000000000000000000e+00,6.666666666666666297e-01,1],[1.000000000000000000e+00,3.333333333333333148e-01,1],[0,0,1],[0,1,1],[1,0,1],[1,1,1],[0.5,1,1],[0.5,0,1],[1.5,0.5,1],[-0.5,0.5,1],[0.5,0.5,1]])
connections = np.array([[0,5],[5,2],[2,6],[6,3],[3,4],[4,1],[1,7],[7,0],[1,8],[8,2],[3,8],[8,0],[4,8],[8,5],[7,8],[8,6],[9,1],[4,9],[8,9],[4,10],[3,10],[8,10],[0,11],[5,11],[8,11],[2,12],[5,12],[8,12],[1,13],[8,13],[7,13],[7,14],[8,14],[0,14],[3,15],[6,15],[8,15],[8,16],[6,16],[2,16],[0,17],[1,18],[2,19],[3,20],[4,21],[5,22],[6,23],[7,24],[8,25],[24,14],[17,14],[25,14],[25,11],[22,11],[17,11],[25,12],[22,12],[19,12],[19,16],[23,16],[25,16],[25,15],[23,15],[20,15],[25,10],[20,10],[21,10],[21,9],[25,9],[18,9],[18,13],[25,13],[24,13],[11,12],[12,16],[16,15],[15,10],[10,9],[9,13],[13,14],[14,11]])

##INOUT SIZE AND SCALE
n_x = 1
n_y = 1
n_z = 2
scale = 1

##NODES AND CONNECTIONS
pattern_nodes = np.zeros((n_x*n_y*n_z*len(nodes), 3))

for i in range(n_x):
    for j in range(n_y):
        for k in range(n_z):
            start_idx = (i*n_y*n_z + j*n_z + k)*len(nodes)
            
            x = nodes[:,0]*scale + i*scale
            y = nodes[:,1]*scale + j*scale
            z = nodes[:,2]*scale + k*scale

            if n_z % 2 == 0:
                if k % 2 == 1:
                    z = (n_z-(k+1))*scale - nodes[:,2]*scale
            else:
                if k % 2 == 1:
                    z = (n_z-(k))*scale - nodes[:,2]*scale

            pattern_nodes[start_idx:start_idx+len(nodes),0] = x
            pattern_nodes[start_idx:start_idx+len(nodes),1] = y
            pattern_nodes[start_idx:start_idx+len(nodes),2] = z
            
pattern_connections = np.zeros((n_x*n_y*n_z*len(connections) + 2*(n_x-1)*(n_y), 2), dtype = int)

unique_nodes = set()
for i in range(len(pattern_nodes)):
    unique_nodes.add(tuple(pattern_nodes[i]))

unique_nodes = np.array(list(unique_nodes))

for i in range(n_x):
    for j in range (n_y):
        for k in range(n_z):
            start_idx = (i*n_y*n_z + j*n_z + k)*len(nodes)

            if j < n_y - 1:
                pattern_connections[2*(i-1)*(n_y) + j, 0] = start_idx + 9
                pattern_connections[2*(i-1)*(n_y) + j, 1] = start_idx + len(nodes) + 11

                pattern_connections[2*(i-1)*(n_y) + j + n_y, 0] = start_idx + 10
                pattern_connections[2*(i-1)*(n_y) + j + n_y, 1] = start_idx + len(nodes) + 12

            for l in range(len(connections)):
                pattern_connections[(i*n_y*n_z + j*n_z + k)*len(connections)+l,0] = start_idx + connections[l,0]
                pattern_connections[(i*n_y*n_z + j*n_z + k)*len(connections)+l,1] = start_idx + connections[l,1]

for i in range(n_x):
    for j in range (n_y):
        for k in range(n_z):
            start_idx = (i*n_y*n_z + j*n_z + k)*len(nodes)
            for connection in connections:
                new_connection = connection + start_idx
                pattern_connections = np.vstack((pattern_connections, new_connection))

##PLOT 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(pattern_nodes[:,0], pattern_nodes[:,1], pattern_nodes[:,2])
# Remove gridlines
ax.grid(False)
# Set white background color
ax.set_facecolor('white')
fig.patch.set_facecolor('white')
ax.set_xlabel('LENGTH')
ax.set_ylabel('WIDTH')
ax.set_zlabel('HEIGHT')
for connections in pattern_connections:
    plt.plot(pattern_nodes[connections,0],pattern_nodes[connections,1],pattern_nodes[connections,2],'k-',lw=0.5)

##TOP VIEW
fig1 = plt.figure()
ax2 = fig1.add_subplot(111)
ax2.scatter(pattern_nodes[:, 0], pattern_nodes[:, 1], c='r')
for connections in pattern_connections:
    plt.plot(pattern_nodes[connections,0],pattern_nodes[connections,1], 'k-')
#fig1.savefig("TOP3DPATTERN.png", dpi=600)

##FRONT VIEW
fig2 = plt.figure()
ax3 = fig2.add_subplot(111)
ax3.scatter(pattern_nodes[:, 1], pattern_nodes[:, 2], c='r')
for connections in pattern_connections:
    plt.plot(pattern_nodes[connections,1],pattern_nodes[connections,2], 'k-')
#fig2.savefig("FRONTOP3DPATTERN.png", dpi=600)

##SIDE VIEW
fig3 = plt.figure()
ax4 = fig3.add_subplot(111)
ax4.scatter(pattern_nodes[:, 0], pattern_nodes[:, 2], c='r')
for connections in pattern_connections:
    plt.plot(pattern_nodes[connections,0],pattern_nodes[connections,2], 'k-')
#fig3.savefig("SIDE3DPATTERN.png", dpi=600)

##SAVE COORIDNATES AND CONNECTION IN LIST
np.savetxt("pattern_nodes.txt", unique_nodes, delimiter=",")
plt.show()

