#4D BODY FITTING MESH CODE FOR JUST A UNIFROM DISTRIBUTED FORCE. #SOME STUFF IS IN ITALIAN SORRY

#POISSONS EQUATION
import numpy as np
import time
import pyamg

start = time.time()

# Constants for the trampoline problem
m = 10
g = 9.81
k = 0.1
forcing_term_constant = m * g / k

def poisson_solver(f, dx, dy, dz, dw, nx, ny, nz, nw, num_levels):
    """
    Solve Poisson's equation on a 4D grid using the finite difference method with cell-centered scheme and multigrid.
    Assumes a zero gradient at the boundaries.
    """

    # Create the cell-centered grid
    x = np.linspace(dx/2, 2 - dx/2, nx)
    y = np.linspace(dy/2, 2.309401 - dy/2, ny)
    z = np.linspace(dz/2, 1.632993 - dz/2, nz)
    w = np.linspace(dw/2, 1.632993 - dw/2, nw)

    A = pyamg.gallery.poisson((nx, ny, nz, nw), format='csr')

    # Create the multigrid solver
    ml = pyamg.smoothed_aggregation_solver(A, max_coarse=num_levels)

    # Solve the system
    b = dx * dy * dz * dw * f.flatten()
    phi = ml.solve(b, tol=1e-10, maxiter=200)

    # Reshape the solution
    phi = phi.reshape((nx, ny, nz, nw))

    phi[0, :, :, :] = 0  # lower x edge
    phi[-1, :, :, :] = 0  # upper x edge
    phi[:, 0, :, :] = 0  # lower y edge
    phi[:, -1, :, :] = 0  # upper y edge
    phi[:, :, 0, :] = 0
    phi[:, :, -1, :] = 0
    phi[:, :, :, 0] = 0
    phi[:, :, :, -1] = 0

    return x, y, z, w, phi

# Define the grid size
nx, ny, nz, nw = 70, 70, 20, 20
# Define the step size of the grid
dx, dy, dz, dw = 2 / (nx - 1), 2.309401 / (ny - 1), 1.632993 / (nz - 1), 1.632993 / (nw - 1)

# Define the function f and change magnitude
f = np.ones((nx, ny, nz, nw))* forcing_term_constant
f[0, :, :, :] = 0
f[-1, :, :, :] = 0
f[:, 0, :, :] = 0
f[:, -1, :, :] = 0
f[:, :, 0, :] = 0
f[:, :, -1, :] = 0
f[:, :, :, 0] = 0
f[:, :, :, -1] = 0


# Solve Poisson's equation with multigrid
num_levels = 16
#x, y, z, w, phi  = poisson_solver(f, dx, dy, dz, dw, nx, ny, nz, nw, num_levels)

#I SCALED THE PROBLEM IN LINE 92 SO THE CONTINUE OF POISSONS EQUATION IS IN LINE 92 TO 116

#IMPORTING 4D MESH FILE. THIS CURRENTLY IS MESHSSK3 WHICH IS THE ONE MODIFIED BY ME.
import plotly.graph_objs as go
import pyvista as pv
import numpy as np
import math
from numpy.linalg import det
import itertools

#this part exports the mesh file and makes a dictionary to put simplices and node coordinates. each simplex and node has an ID
# Leggi i dati dal file .txt
##CHANGE INPUT FILE NAME TO ORIGINAL4DMESH.txt IF YOU WANT TO USE THE ORIGINAL MESH WITH DISCONNECTED NODES
with open('NEW4DMESH.txt', 'r') as f:
    lines = f.readlines()
# Salva i dati in due liste separate
## FOR ORIGINAL MESH
# simplices = [list(map(int, l.strip().split())) for l in lines[5:3557]]
# nodes = [list(map(float, l.strip().split())) for l in lines[3561:4032]]
##FOR NEW MESH
simplices = [list(map(int, l.strip().split())) for l in lines[5:3701]]
nodes = [list(map(float, l.strip().split())) for l in lines[3705:4176]]
# Crea un dizionario per il mesh 4D

# Minimum and maximum coordinates of the 4D mesh nodes
min_coords = np.min(nodes, axis=0)[1:]  # Ignore the first column, which contains node IDs
max_coords = np.max(nodes, axis=0)[1:]

desired_range = (1, 1000)
scaling_factor = (desired_range[1] - desired_range[0]) / (max_coords - min_coords)

# Scaling step sizes (dx, dy, dz, dw)
scaled_dx = dx * scaling_factor[0]
scaled_dy = dy * scaling_factor[1]
scaled_dz = dz * scaling_factor[2]
scaled_dw = dw * scaling_factor[3]

f = np.ones((nx, ny, nz, nw)) * forcing_term_constant
f[0, :, :, :] = 0
f[-1, :, :, :] = 0
f[:, 0, :, :] = 0
f[:, -1, :, :] = 0
f[:, :, 0, :] = 0
f[:, :, -1, :] = 0
f[:, :, :, 0] = 0
f[:, :, :, -1] = 0

# Solve Poisson's equation with multigrid and scaled step sizes
x, y, z, w, phi = poisson_solver(f, scaled_dx, scaled_dy, scaled_dz, scaled_dw, nx, ny, nz, nw, num_levels)

mesh = {}
# Crea una lista di nodi vuota
nodes_list = [None] * len(nodes)
# Aggiungi i nodi al dizionario e alla lista di nodi
min_coords = np.min(nodes, axis=0)[1:]  # Ignore the first column, which contains node IDs
max_coords = np.max(nodes, axis=0)[1:]

for node in nodes:
    node_id = int(node[0])
    coords = tuple((np.array(node[1:]) - min_coords) * scaling_factor + desired_range[0])
    mesh[node_id] = {'coords': coords}
    nodes_list[node_id - 1] = [node_id] + list(coords)
# for node in nodes:
#     node_id = int(node[0])
#     coords = tuple(node[1:])
#     mesh[node_id] = {'coords': coords}
#     nodes_list[node_id - 1] = node  #lista di tutti i nodi con indici file

# Aggiungi i simplessi al dizionario
for simplex in simplices:
    simplex_id = simplex[0]
    node_ids = simplex[1:]
    tet_coords = tuple(mesh[node_id]['coords'] for node_id in node_ids)
    mesh[simplex_id] = {'type': 'tetrahedron', 'nodes': node_ids, 'coords': tet_coords}

#EXAMPLE OF HOW TO ACCESS THE NODE COORIDNATES, OF THE SIMPLEX, WITH INPUT. JUST REMOVER THE "#" FROM LINE 150 TO 153 TO PRINT IT.

# Esempio di come accedere alle coordinate di un nodo e le coordinate dei nodi di un simplex
simplex_id = 10
simplex_node_ids = mesh[simplex_id]['nodes']
simplex_node_coords = [nodes_list[node_id - 1] for node_id in simplex_node_ids]
# print(f"I nodi del simplex con id {simplex_id} sono: {simplex_node_ids}")
# print(f"Le coordinate dei nodi del simplex con id {simplex_id} sono:")
# for node in simplex_node_coords:
#     print(node)

#Example to how access the simplices around a node
def simplices_around_node(node_id, mesh):
    connected_simplices = []
    for simplex_id, simplex_info in mesh.items():
        if simplex_info['type'] == 'tetrahedron' and node_id in simplex_info['nodes']:
            connected_simplices.append(simplex_id)
    return connected_simplices
node_id_to_print = 2
connected_simplices = simplices_around_node(node_id_to_print, mesh)
#print(f"Simplices around node {node_id_to_print}: {connected_simplices}")

#print(mesh)

from collections import defaultdict
connected_simplices = defaultdict(list)
for simplex in simplices:
    for i in simplex[1:]:
        connected_simplices[i].append(simplex)

#VOLUMES of simplices calculation
volumes = []
for simplex in simplices:
    simplex_id = simplex[0]
    node_ids = simplex[1:]
    simplex_node_coords = [nodes_list[node_id - 1][1:] for node_id in node_ids]  # extract the node coordinates
    mat = np.hstack((simplex_node_coords, np.ones((5,1))))  # create the matrix
    volume = 1/np.math.factorial(4) * np.abs(np.linalg.det(mat))  # compute the volume using the formula
    volumes.append(volume)

# compute centroids of simplices
centroids = {}
for simplex in simplices:
    simplex_id = simplex[0]
    node_ids = simplex[1:]
    simplex_node_coords = [nodes_list[node_id - 1][1:] for node_id in node_ids] 
    centroid_coords = np.mean(simplex_node_coords, axis=0)
    centroids[simplex_id] = {'coords': centroid_coords}

# compute medians of edges of simplices
edge_medians = {}
for simplex in simplices:
    simplex_id = simplex[0]
    node_ids = simplex[1:]
    for i, j in itertools.combinations(node_ids, 2):
        edge_id = (i, j) if i < j else (j, i)
        if edge_id not in edge_medians:
            node1_coords = np.array(nodes_list[i - 1][1:])
            node2_coords = np.array(nodes_list[j - 1][1:])
            median_coords = 0.5 * (node1_coords + node2_coords)
            edge_medians[edge_id] = {'coords': median_coords, 'simplices': []}
        edge_medians[edge_id]['simplices'].append(simplex_id)

# Get the node ids of simplex, example of accessing the edge emdians of a simplex
simplex_id = 2
node_ids = mesh[simplex_id]['nodes']
# Get the coordinates of the edge medians of simplex 1
edge_median_coords = []
for i, j in itertools.combinations(node_ids, 2):
    edge_id = (i, j) if i < j else (j, i)
    edge_median_coords.append(edge_medians[edge_id]['coords'])
# # Print the coordinates of the edge medians of simplex
# print("Edge median coordinates of simplex 1:")
# for coord in edge_median_coords:
#     print(coord)

# compute centroids of faces of each simplex
face_centroids = {}
for simplex in simplices:
    simplex_id = simplex[0]
    node_ids = simplex[1:]
    for i, j, k in itertools.combinations(node_ids, 3):
        face_id = tuple(sorted((i, j, k)))
        if face_id not in face_centroids:
            node1_coords = nodes_list[i - 1][1:]
            node2_coords = nodes_list[j - 1][1:]
            node3_coords = nodes_list[k - 1][1:]
            centroid_coords = np.mean([node1_coords, node2_coords, node3_coords], axis=0)
            face_centroids[face_id] = {'coords': centroid_coords, 'simplices': []}
        face_centroids[face_id]['simplices'].append(simplex_id)

#Function to calculate volume of a simplex determinant matrix
def volume_4d_simplex(simplex_coords):
    # Get the 5 vertices of the simplex
    p1, p2, p3, p4, p5 = simplex_coords
    
    # Create the matrix
    mat = np.array([
        [p1[0], p1[1], p1[2], p1[3], 1],
        [p2[0], p2[1], p2[2], p2[3], 1],
        [p3[0], p3[1], p3[2], p3[3], 1],
        [p4[0], p4[1], p4[2], p4[3], 1],
        [p5[0], p5[1], p5[2], p5[3], 1]
    ])
    
    # Compute the volume using the formula
    volume = 1/np.math.factorial(4) * np.abs(np.linalg.det(mat))
    return volume

# Calculate the volume of the median dual cells (CONTROL VOLUMES OF EACH NODE)
median_dual_cell_volumes = np.zeros(len(nodes_list))
#median_dual_cell_volumes = {}

for simplex in simplices:
    simplex_id = simplex[0]
    node_ids = simplex[1:]
    
    # Get the centroid of the current simplex
    centroid_coords = centroids[simplex_id]['coords']
    
    for i in node_ids:
        node_coords = nodes_list[i - 1][1:]
        
        # Get the centroids of the faces connected to the current node
        connected_face_centroids = [face_centroids[face_id]['coords'] for face_id in face_centroids if i in face_id and simplex_id in face_centroids[face_id]['simplices']]
        
        # Get the medians of the edges connected to the current node
        connected_edge_medians = [edge_medians[edge_id]['coords'] for edge_id in edge_medians if i in edge_id and simplex_id in edge_medians[edge_id]['simplices']]
        
        # Calculate the volume of the median dual cell
        median_dual_cell_volume = 0
        for face_centroid in connected_face_centroids:
            for edge_median1, edge_median2 in itertools.combinations(connected_edge_medians, 2):
                simplex_coords = (node_coords, edge_median1, edge_median2, face_centroid, centroid_coords)
                median_dual_cell_volume += volume_4d_simplex(simplex_coords)
        
        median_dual_cell_volumes[i - 1] += median_dual_cell_volume
        #median_dual_cell_volumes[node_id] = median_dual_cell_volumes.get(node_id, 0) + median_dual_cell_volume

#TO PRINT MEDIAN DUAL CELL VOLUMES
#print("Median dual cell volumes:", median_dual_cell_volumes)


#DECOMPOSITION OF SIMPLEX IN POLYTOPES TO FIND COMMON FACES BETWEEN THEM AND SO BETWEEN TWO NODE
def decompose_simplex_into_polytopes(simplex_id, mesh, face_centroids, edge_medians, centroids):
    node_ids = mesh[simplex_id]['nodes']
    simplex_centroid = centroids[simplex_id]['coords']
    polytopes = []
    for node_id in node_ids:
        node_coords = mesh[node_id]['coords']
        
        face_ids = [face_id for face_id in face_centroids if node_id in face_id and simplex_id in face_centroids[face_id]['simplices']]
        face_centroids_coords = [face_centroids[face_id]['coords'] for face_id in face_ids]
        
        edge_ids = [edge_id for edge_id in edge_medians if node_id in edge_id and simplex_id in edge_medians[edge_id]['simplices']]
        edge_medians_coords = [edge_medians[edge_id]['coords'] for edge_id in edge_ids]
        
        polytope_vertices = [node_coords] + face_centroids_coords + edge_medians_coords + [simplex_centroid]
        polytopes.append(polytope_vertices)
    
    return polytopes

def find_common_faces_between_polytopes(polytopes):
    common_faces = []
    for polytope1, polytope2 in itertools.combinations(polytopes, 2):
        intersection = set(tuple(sorted(vertex, key=lambda x: str(x))) for vertex in polytope1).intersection(
            set(tuple(sorted(vertex, key=lambda x: str(x))) for vertex in polytope2))
        if len(intersection) >= 5:  # a common 3D face should have at least 4 vertices
            common_faces.append(list(intersection))
    return common_faces

node_info = {node_id: {'faces': [], 'cv_face_centroids': [], 'neighbors': [], 'face_areas': []} for node_id in range(1, len(nodes_list) + 1)}
for simplex in simplices:
    simplex_id = simplex[0]
    polytopes = decompose_simplex_into_polytopes(simplex_id, mesh, face_centroids, edge_medians, centroids)
    common_faces = find_common_faces_between_polytopes(polytopes)

    # Update the node_info dictionary here
    for i, polytope1 in enumerate(polytopes):
        node1_id = mesh[simplex_id]['nodes'][i]
        for j, polytope2 in enumerate(polytopes):
            node2_id = mesh[simplex_id]['nodes'][j]

            # Skip if the polytopes belong to the same node
            if node1_id == node2_id:
                continue

            # Find the common face between polytope1 and polytope2
            common_face = set(map(tuple, polytope1)).intersection(set(map(tuple, polytope2)))
            if len(common_face) >= 5:  # a common 4D face should have at least 5 vertices
                common_face = [np.array(vertex) for vertex in common_face]

                # Calculate the centroid of the common face
                cv_face_centroid = np.mean(common_face, axis=0)
                if len(common_face) == 5:
                    face_area = volume_4d_simplex(common_face)
                else:
                    print("Unexpected number of vertices:", len(common_face))
                    print("common_face:", common_face)

                    # Update the node_info dictionary for node1_id and node2_id
                if node2_id not in node_info[node1_id]['neighbors']:
                    node_info[node1_id]['faces'].append(common_face)
                    node_info[node1_id]['cv_face_centroids'].append(cv_face_centroid)
                    node_info[node1_id]['neighbors'].append(node2_id)
                    node_info[node1_id]['face_areas'].append(face_area)

                if node1_id not in node_info[node2_id]['neighbors']:
                    node_info[node2_id]['faces'].append(common_face)
                    node_info[node2_id]['cv_face_centroids'].append(cv_face_centroid)
                    node_info[node2_id]['neighbors'].append(node1_id)
                    node_info[node2_id]['face_areas'].append(face_area)
node_id = 200
print(f"Node {node_id} has the following information:")
print(f"  Faces: {node_info[node_id]['faces']}")
print(f"  cv_Face centroids: {node_info[node_id]['cv_face_centroids']}")
print(f"  Neighboring nodes: {node_info[node_id]['neighbors']}")
print(f"  Faces Areas: {node_info[node_id]['face_areas']}")

def calculate_distances_to_face_centroids(node_info, nodes_list):
    for node_id, node_data in node_info.items():
        node_coords = np.array(nodes_list[node_id - 1][1:])
        distances = []
        neighbor_distances = []
        for i, cv_face_centroid in enumerate(node_data['cv_face_centroids']):
            distance = np.linalg.norm(node_coords - np.array(cv_face_centroid))
            distances.append(distance)
            neighbor_id = node_data['neighbors'][i]
            neighbor_coords = np.array(nodes_list[neighbor_id - 1][1:])
            neighbor_distance = np.linalg.norm(neighbor_coords - np.array(cv_face_centroid))
            neighbor_distances.append(neighbor_distance)
        node_info[node_id]['distances'] = distances
        node_info[node_id]['neighbor_distances'] = neighbor_distances

calculate_distances_to_face_centroids(node_info, nodes_list)

#print(f"Node {node_id} distances to face centroids: {node_info[node_id]['distances']}")
#print(f"Node {node_id} neighboring nodes distances to face centroids: {node_info[node_id]['neighbor_distances']}")

#selecting mass and constant for each node
mass = 1 # kg
gamma = 1.4

# Update nodes_list with volumes and masses and calculate pressure due to gas law in each node
for i, node in enumerate(nodes_list):
    volume = median_dual_cell_volumes[i]
    density = mass / volume
    pressure_grad = (density ** gamma)
    node.extend([volume, density])
    node.append(pressure_grad)

#storing the pressure grad, volume and mass in node id information
#node_id = 470
volume = nodes_list[node_id - 1][5]
density = nodes_list[node_id - 1][6]
pressure = nodes_list[node_id - 1][7]


#INTERPOLATION OF POISSON S EQUATION TO 4D GRID
def interpolate_phi_at_centroid(phi, x, y, z, w, c_centroid):
    """
    Interpolate the value of phi at the given 4D mesh node coordinates using tetra-linear interpolation.

    :param phi: The Poisson's solver solution (phi).
    :param x: The x-coordinates of the cell-centered grid.
    :param y: The y-coordinates of the cell-centered grid.
    :param z: The z-coordinates of the cell-centered grid.
    :param w: The w-coordinates of the cell-centered grid.
    :param coords: The 4D mesh node coordinates.
    :return: The interpolated value of phi at the given 4D mesh node coordinates.
    """
    xi, yi, zi, wi = c_centroid
    idx_x = np.searchsorted(x, xi) - 1
    idx_y = np.searchsorted(y, yi) - 1
    idx_z = np.searchsorted(z, zi) - 1
    idx_w = np.searchsorted(w, wi) - 1

    dx = x[1] - x[0]
    dy = y[1] - y[0]
    dz = z[1] - z[0]
    dw = w[1] - w[0]

    tx = (xi - x[idx_x]) / dx
    ty = (yi - y[idx_y]) / dy
    tz = (zi - z[idx_z]) / dz
    tw = (wi - w[idx_w]) / dw

    interpolated_value = 0
    for i, j, k, l in itertools.product(range(2), repeat=4):
        x_idx = idx_x + i
        y_idx = idx_y + j
        z_idx = idx_z + k
        w_idx = idx_w + l

        # Ensure the indices are within the valid range
        if (0 <= x_idx < phi.shape[0]) and (0 <= y_idx < phi.shape[1]) and \
           (0 <= z_idx < phi.shape[2]) and (0 <= w_idx < phi.shape[3]):

            weight = ((1 - tx) if i == 0 else tx) * \
                     ((1 - ty) if j == 0 else ty) * \
                     ((1 - tz) if k == 0 else tz) * \
                     ((1 - tw) if l == 0 else tw)
            interpolated_value += weight * phi[x_idx, y_idx, z_idx, w_idx]

    #interpolated_value = np.clip(interpolated_value, phi.min(), phi.max())

    return interpolated_value

#STORE INTERPOLATED VALUE IN CV FACE CENTROIDS
def store_phi_at_cv_face_centroids(node_info, phi, x, y, z, w):
    for node_id, node_data in node_info.items():
        cv_face_centroids = node_data['cv_face_centroids']
        centroid_phi_values = []
        for c_centroid in cv_face_centroids:
            centroid_phi = interpolate_phi_at_centroid(phi, x, y, z, w, c_centroid)
            centroid_phi_values.append(centroid_phi)
        node_info[node_id]['centroid_phi_values'] = centroid_phi_values

store_phi_at_cv_face_centroids(node_info, phi, x, y, z, w)


#GEOMETRIC INTERPOLATION OF GAS LAW PRESSURE TO FACE CENTROIDS USING PHD STUDENT EQUATION
def calculate_interface_pressures_and_store(node_id, node_info, nodes_list):
    # Retrieve node data
    distances = node_info[node_id]['distances']
    neighbor_distances = node_info[node_id]['neighbor_distances']
    neighbors = node_info[node_id]['neighbors']
    pressure = nodes_list[node_id - 1][7]
    centroid_phi_values = node_info[node_id]['centroid_phi_values']

    # Calculate the geometric interpolation factor and interface pressures
    interface_pressures = []
    for i, neighbor_id in enumerate(neighbors):
        neighbor_pressure = nodes_list[neighbor_id - 1][7]
        dCf = distances[i]
        dfF = neighbor_distances[i]

        # Calculate the geometric interpolation factor (gf)
        gf = dCf / (dCf + dfF)

        # Calculate the interface pressure
        Pgi_half = gf * neighbor_pressure + (1 - gf) * pressure
        Pgi_half += centroid_phi_values[i]  # Add phi value from Poisson's equation
        interface_pressures.append(Pgi_half)

    # Store interface pressures in node_info
    node_info[node_id]['interface_pressures'] = interface_pressures

calculate_interface_pressures_and_store(node_id, node_info, nodes_list)


#FORCES, THIS PART IS A BIT EXPERIMENTAL 
def calculate_forces(node_info):
    for node_id, node_data in node_info.items():
        interface_pressures = node_data['interface_pressures']
        face_areas = node_data['face_areas']
        forces = []

        for i, (pressure, area) in enumerate(zip(interface_pressures, face_areas)):
            force = -pressure * area
            forces.append(force)

        node_info[node_id]['forces'] = forces

for node_id in range(1, len(nodes_list) + 1):
    calculate_interface_pressures_and_store(node_id, node_info, nodes_list)

calculate_forces(node_info)
node_id = 200
print(f"Interface pressures for node {node_id}: {node_info[node_id]['interface_pressures']}")
print(f"Interface forces for node {node_id}: {node_info[node_id]['forces']}")

print(f'duration: {time.time() - start}')

# Calculate speed of sound and local time step for each node
local_time_steps = []
for node_id, node_data in node_info.items():
    node_volume = nodes_list[node_id - 1][5]
    node_pressure = nodes_list[node_id - 1][7]
    node_density = nodes_list[node_id - 1][6]
    node_velocity = node_volume / node_density

    # Calculate displacement limits
    distances = node_info[node_id]['distances']

    if len(distances) > 0:
        delta_x = min(distances)
    else:
    # Handle the case when there are no distances available
    # You can skip this node or handle it in another way
        print(f"No distances found for node {node_id}")
        continue  # skips to the next iteration of the loop

    # Calculate speed of sound
    speed_of_sound = np.sqrt(gamma * node_pressure / node_density)

    # Calculate local time step
    #dt_local = delta_x / (GAMMA * node_pressure * node_velocity + speed_of_sound)
    dt_local = delta_x / abs(node_velocity + speed_of_sound)
    local_time_steps.append(dt_local)

# Calculate global time step
dt_global = min(local_time_steps)

# Update node velocities and positions
for node_id, node_data in node_info.items():
    node_forces = np.array(node_data['forces'])
    node_mass = nodes_list[node_id - 1][6] * nodes_list[node_id - 1][5]

    # Calculate acceleration and velocity of the node
    acceleration = node_forces.sum(axis=0) / node_mass
    velocity = acceleration * dt_global

    # Update node position
    node_coords = np.array(nodes_list[node_id - 1][1:5])
    new_node_coords = node_coords + velocity * dt_global

    # Update the nodes_list with the new position
    nodes_list[node_id - 1][1:5] = new_node_coords.tolist()

#SOME VISUALISATION TRIALS
# # Assuming nodes_list contains the coordinates of all the nodes
# mesh_data = np.array([node_info[1:5] for node_info in nodes_list])
# from sklearn.decomposition import PCA

# # Initialize PCA object to reduce the dimensions to 3
# pca = PCA(n_components=3)

# # Fit and transform the mesh_data using PCA
# reduced_mesh_data = pca.fit_transform(mesh_data)

# # Visualize the reduced data using Matplotlib
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# ax.scatter(reduced_mesh_data[:, 0], reduced_mesh_data[:, 1], reduced_mesh_data[:, 2])

# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')

# #plt.show()

# #print(new_node_coords)

# # Write the updated node coordinates and IDs to a text file
# with open('updated_nodes.txt', 'w') as output_file:
#     for node in nodes_list:
#         node_id = node[0]
#         coordinates = node[1:5]
#         output_file.write(f"{node_id} {coordinates[0]} {coordinates[1]} {coordinates[2]} {coordinates[3]}\n")


# @article{BeOlSc2022,
#   author    = {Nathan Bell and Luke N. Olson and Jacob Schroder},
#   title     = {{PyAMG}: Algebraic Multigrid Solvers in Python},
#   journal   = {Journal of Open Source Software},
#   year      = {2022},
#   publisher = {The Open Journal},
#   volume    = {7},
#   number    = {72},
#   pages     = {4142},
#   doi       = {10.21105/joss.04142},
#   url       = {https://doi.org/10.21105/joss.04142},
# }



