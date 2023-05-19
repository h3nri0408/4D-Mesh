##2D PATTERN CREATION

import matplotlib.pyplot as plt
import numpy as np

nodes = np.array([[0,0],[1,0],[2,0],[3,0],[4,0],[0.5,0.2],[1.5,0.2],[2.5,0.2],[3.5,0.2],[0,0.2],[1,0.2],[2,0.2],[3,0.2],[4,0.2],[0,0.4],[1,0.4],[2,0.4],[3,0.4],[4,0.4]])

connections = np.array([[0,1],[1,2],[2,3],[3,4],[0,5],[5,1],[1,6],[6,2],[2,7],[7,3],[3,8],[8,4],[0,9],[1,10],[2,11],[3,12],[4,13],[9,10],[10,11],[11,12],[12,13],[9,14],[10,15],[11,16],[12,17],[13,18],[14,15],[15,16],[16,17],[17,18],[14,5],[5,15],[15,6],[6,16],[16,7],[7,17],[17,8],[8,18],[5,10],[6,10]])

##INPUT SIZE AND SCALE
n_x = 3
n_y = 3
scale = 1

##NODES AND CONNECTIONS
pattern_nodes = np.zeros((n_x * n_y * len(nodes), 2))

for i in range(n_x):
    for j in range(n_y):
        start_idx = (i * n_y + j) * len(nodes)

        x = nodes[:, 0] * scale + i * scale
        y = nodes[:, 1] * scale + j * (nodes[:, 1].max() - nodes[:, 1].min())

        pattern_nodes[start_idx:start_idx + len(nodes), 0] = x
        pattern_nodes[start_idx:start_idx + len(nodes), 1] = y

pattern_connections = np.zeros((n_x * n_y * len(connections), 2), dtype=int)

for i in range(n_x):
    for j in range(n_y):
        start_idx = (i * n_y + j) * len(nodes)

        for k in range(len(connections)):
            pattern_connections[(i * n_y + j) * len(connections) + k, 0] = start_idx + connections[k, 0]
            pattern_connections[(i * n_y + j) * len(connections) + k, 1] = start_idx + connections[k, 1]

##PLOT 
fig = plt.figure()
plt.scatter(pattern_nodes[:, 0], pattern_nodes[:, 1], color='red')

#yellow_connections_layer1 = [30, 31, 32, 33, 38, 39, 40, 22]  # Specify the indices of the connections you want to change the color for in the first layer.
yellow_connections_layer1 = []
#yellow_connections_layer2 = [5, 6, 4, 7, 0, 1, 38, 39, 40, 13]  # Specify the indices of the connections you want to change the color for in the second layer.
yellow_connections_layer2 = []

for i in range(n_x):
    for j in range(n_y):
        for k in range(len(connections)):
            start_idx = (i * n_y + j) * len(connections)
            connections = pattern_connections[start_idx:start_idx + len(connections)]

            if j == 0 and k in yellow_connections_layer1:
                color = 'y'
            elif j == 1 and k in yellow_connections_layer2:
                color = 'y'
            else:
                color = 'b'

            plt.plot(pattern_nodes[connections[k], 0], pattern_nodes[connections[k], 1], f'{color}-')

plt.ylim([-0.1, n_y * (nodes[:, 1].max() - nodes[:, 1].min()) + 0.1])
plt.show()
#fig.savefig("LINES6.png", dpi=600)
