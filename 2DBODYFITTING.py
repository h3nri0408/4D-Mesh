import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Delaunay
from collections import defaultdict
from scipy import interpolate
from scipy.interpolate import interp2d

##CREATION OF 2D PATTTERN
nodes = np.array([[0,0],[0,1],[1,0],[1,1],[0.5,1],[0.5,0],[1.5,0.5],[-0.5,0.5],[0.5,0.5]])
connections = np.array([[0,5],[5,2],[2,6],[6,3],[3,4],[4,1],[1,7],[7,0],[1,8],[8,2],[3,8],[8,0],[4,8],[8,5],[7,8],[8,6]])

n_x = 5
n_y = 5
scale = 0.5

# Delaunay triangulation from nodes
tri = Delaunay(nodes)

pattern_nodes = []

for i in range(n_x):
    for j in range(n_y):
        x = nodes[:,0]*scale + i*scale
        y = nodes[:,1]*scale + j*scale
        pattern_nodes.extend(list(zip(x, y)))

pattern_nodes = np.array(pattern_nodes)

pattern_connections = []

for i in range(n_x):
    for j in range(n_y):
        start_idx = (i*n_y + j)*len(nodes)
        for k in range(len(connections)):
            pattern_connections.append((start_idx + connections[k, 0], start_idx + connections[k, 1]))

pattern_connections = np.array(pattern_connections)

unique_nodes, unique_indices = np.unique(pattern_nodes, axis=0, return_inverse=True)

new_pattern_connections = []

for connection in pattern_connections:
    new_pattern_connections.append((unique_indices[connection[0]], unique_indices[connection[1]]))

pattern_connections = np.array(new_pattern_connections)

# Remove duplicate connections
unique_connections = {tuple(sorted(conn)) for conn in pattern_connections}
pattern_connections = np.array(list(unique_connections))
fig = plt.figure()
plt.scatter(unique_nodes[:,0], unique_nodes[:,1])
plt.xlabel('LENGTH')
plt.ylabel('WIDTH')

for connection in pattern_connections:
    plt.plot(unique_nodes[connection,0], unique_nodes[connection,1], 'k-')
fig.savefig("2DPATTERN.png", dpi=600)
#plt.show()

# print("Unique nodes:")
# for i, node in enumerate(unique_nodes):
#     print(f"Node {i}: {node}")

# print("\nPattern connections:")
# for i, connection in enumerate(pattern_connections):
#     print(f"Connection {i}: {connection}")

# Find triangles
triangles = set()
for conn1 in new_pattern_connections:
    for conn2 in new_pattern_connections:
        if not np.array_equal(conn1, conn2):
            # Find the common node between the two connections
            common_nodes = set(conn1).intersection(conn2)
            if len(common_nodes) == 1:
                common_node = common_nodes.pop()

                # Find the other two nodes of the potential triangle
                node1_set = set(conn1).difference({common_node})
                node2_set = set(conn2).difference({common_node})
                if node1_set and node2_set:
                    node1 = node1_set.pop()
                    node2 = node2_set.pop()

                    # Check if there is a connection between the two other nodes
                    if {node1, node2} in [set(conn) for conn in new_pattern_connections]:
                        # Add the triangle to the set, sorting the node indices to avoid duplicates
                        triangles.add(tuple(sorted([common_node, node1, node2])))


## GENERAL CHECKS TO SEE IF EVERYTHING LOOKS OK
# print("\nNumber of triangles:", len(triangles))

# # Create the nodes list with nodes and their indices
# nodes_list = [(idx, tuple(node)) for idx, node in enumerate(unique_nodes)]
# print(len(nodes_list))
# # Create the triangles list with triangle IDs and the 3 nodes they are composed of
# triangles_list = [(idx, [unique_nodes[node_idx] for node_idx in triangle]) for idx, triangle in enumerate(triangles)]
# print(len(triangles_list))
# # Print the nodes list
# print("\nNodes list:")
# for node in nodes_list:
#     print(f"Node index: {node[0]}, Node coordinates: {node[1]}")
# # Print the triangles list
# print("\nTriangles list:")
# for triangle in triangles_list:
#     print(f"Triangle ID: {triangle[0]}, Nodes: {triangle[1]}")
# print(nodes_list)
# print(triangles_list)

# Create the nodes list with nodes and their indices
nodes_list = [[idx] + node.tolist() for idx, node in enumerate(unique_nodes)]

# Create the triangles list with triangle IDs and the 3 node ids they are composed of
triangles_list = [[idx] + list(triangle) for idx, triangle in enumerate(triangles)]

# Calculate the centroids of the triangles
centroids = []
for triangle in triangles_list:
    centroid = np.mean([unique_nodes[triangle[1]], unique_nodes[triangle[2]], unique_nodes[triangle[3]]], axis=0)
    centroids.append(centroid)

# Calculate the medians of the edges of all triangles without repeating the medians for the same edge
median_edges_set = set()
for triangle in triangles_list:
    for i in range(1, 4):
        node1 = unique_nodes[triangle[i]]
        node2 = unique_nodes[triangle[i % 3 + 1]]
        midpoint = tuple((np.array(node1) + np.array(node2)) / 2)

        # Ensure the midpoint is unique by sorting the node IDs
        edge = tuple(sorted((triangle[i], triangle[i % 3 + 1])))

        # Add the edge to the set of median edges
        median_edges_set.add(edge)

# Calculate the midpoints
median_midpoints = {edge: tuple((unique_nodes[edge[0]] + unique_nodes[edge[1]]) / 2) for edge in median_edges_set}
# Plot the centroids
for centroid in centroids:
    plt.plot(centroid[0], centroid[1], 'ro')

# Plot the midpoints
for midpoint in median_midpoints.values():
    plt.plot(midpoint[0], midpoint[1], 'go')
plt.show()




from collections import defaultdict
import numpy as np

# Step 1: Data structure to store the control volume faces for each node
control_volume_faces = defaultdict(list)

# Step 2: Iterate through the triangles and centroids
for triangle, centroid in enumerate(centroids):
    
    # Step 3: For each triangle and centroid, iterate through the nodes of the triangle
    for i in range(1, 4):
        
        # Step 4: For each node, find the corresponding median points and centroid
        node = triangles_list[triangle][i]
        edges = [tuple(sorted((node, triangles_list[triangle][j]))) for j in range(1, 4) if j != i]
        midpoints = [median_midpoints[edge] for edge in edges]
        
        # Step 5: Store the faces as control volume faces for the corresponding node
        control_volume_faces[node].append((centroid, midpoints[0], midpoints[1]))

# Function to calculate the area of a triangle given its vertices
def triangle_area(a, b, c):
    return 0.5 * abs((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]))

# Step 6: Calculate the control volume areas for each node
control_volume_areas = {}
for node_id, faces in control_volume_faces.items():
    area = 0
    for face in faces:
        centroid, midpoint1, midpoint2 = face
        # Split the polygon into two triangles and sum their areas
        area += triangle_area(unique_nodes[node_id], centroid, midpoint1)
        area += triangle_area(unique_nodes[node_id], centroid, midpoint2)
    control_volume_areas[node_id] = area

##CHECKS
# # Print the control volume areas for each node
# for node_id, area in control_volume_areas.items():
#     print(f"Node {node_id}: Control Volume Area = {area}")

# print("Total control volume area:", sum(control_volume_areas.values()))


chosen_node_id = 71  # Change this to visualiSe the control volume for a different node
##PLOT OF CONTROL VOLUME AROUND A NODE IN YELLOW
# Plot the 2D pattern
plt.scatter(unique_nodes[:, 0], unique_nodes[:, 1], color='black', s=5)
plt.scatter(np.array(centroids)[:, 0], np.array(centroids)[:, 1], color='black', s=20)

for connection in pattern_connections:
    plt.plot(unique_nodes[connection, 0], unique_nodes[connection, 1], 'k-')

# Plot connections between centroids and corresponding median edges
for triangle, centroid in enumerate(centroids):
    for i in range(1, 4):
        edge = tuple(sorted((triangles_list[triangle][i], triangles_list[triangle][i % 3 + 1])))
        midpoint = median_midpoints[edge]
        plt.plot([centroid[0], midpoint[0]], [centroid[1], midpoint[1]], 'r-')

# Plot the control volume area around the chosen node
if chosen_node_id in control_volume_faces:
    chosen_node = unique_nodes[chosen_node_id]
    for centroid, midpoint1, midpoint2 in control_volume_faces[chosen_node_id]:
        plt.fill([chosen_node[0], centroid[0], midpoint1[0]], [chosen_node[1], centroid[1], midpoint1[1]], color='yellow', alpha=0.5)
        plt.fill([chosen_node[0], centroid[0], midpoint2[0]], [chosen_node[1], centroid[1], midpoint2[1]], color='yellow', alpha=0.5)
        plt.plot([chosen_node[0], centroid[0]], [chosen_node[1], centroid[1]], 'g--', linewidth=0.5)

plt.show()


# Create a data structure to store the control volume details for each node
node_info = {}

for node_id, faces in control_volume_faces.items():
    # Data structure for this node
    node_info[node_id] = {'CONTROL VOLUME AREA': control_volume_areas[node_id], 
                          'FACES': {}}

    for face in faces:
        centroid, midpoint1, midpoint2 = face

        # Calculate the length of the edges
        length1 = np.linalg.norm(np.array(midpoint1) - np.array(centroid))
        length2 = np.linalg.norm(np.array(midpoint2) - np.array(centroid))

        # Calculate the midpoints of the edges
        midpoint_face1 = tuple(0.5 * (np.array(midpoint1) + np.array(centroid)))
        midpoint_face2 = tuple(0.5 * (np.array(midpoint2) + np.array(centroid)))

        # Convert numpy arrays to tuple
        centroid = tuple(centroid)
        midpoint1 = tuple(midpoint1)
        midpoint2 = tuple(midpoint2)

        # Store the face details
        face_id1 = frozenset([centroid, midpoint1])
        face_id2 = frozenset([centroid, midpoint2])
        node_info[node_id]['FACES'][face_id1] = {'MIDPOINT': midpoint_face1, 'LENGTH': length1}
        node_info[node_id]['FACES'][face_id2] = {'MIDPOINT': midpoint_face2, 'LENGTH': length2}

# Identify the neighbouring nodes and shared faces
for node_id, info in node_info.items():
    for face_id, face_info in info['FACES'].items():
        for other_node_id, other_info in node_info.items():
            if node_id != other_node_id and face_id in other_info['FACES']:
                face_info['NEIGHBOURING NODE'] = other_node_id
##CHECK
# # Print the control volume details for each node
# for node_id, info in node_info.items():
#     print(f"Node {node_id}:")
#     print(f"  Control Volume Area = {info['CONTROL VOLUME AREA']}")
#     for face_id, face_info in info['FACES'].items():
#         print(f"  Face {face_id}:")
#         print(f"    Midpoint = {face_info['MIDPOINT']}")
#         print(f"    Length = {face_info['LENGTH']}")
#         if 'NEIGHBOURING NODE' in face_info:
#             print(f"    Neighbouring node = {face_info['NEIGHBOURING NODE']}")

# Assign a mass to each node
mass = 0.1

# Set the heat capacity ratio
gamma = 1.4

# Calculate and store the pressure for each node
for node_id, info in node_info.items():
    volume = info['CONTROL VOLUME AREA']
    pressure = (mass / volume) ** gamma
    node_info[node_id]['PRESSURE'] = pressure

# # Print the pressure for each node
# for node_id, pressure in node_info.items():
#     print(f"Node {node_id}: Pressure = {pressure}")

# Create arrays of the node coordinates and pressures
node_coords = unique_nodes
node_pressures = np.array([info['PRESSURE'] for info in node_info.values()])

# Create an interpolator for the node pressures
pressure_interpolator = interp2d(node_coords[:,0], node_coords[:,1], node_pressures, kind='linear')


# Interpolate the pressure values at the midpoints of the faces
for node_id, info in node_info.items():
    for face_id, face_info in info['FACES'].items():
        # Calculate the interpolated pressure value
        interpolated_pressure = pressure_interpolator(*face_info['MIDPOINT'])
        # Store the interpolated pressure value
        face_info['PRESSURE'] = interpolated_pressure
##CHECK
# #print the details for each node
# for node_id, info in node_info.items():
#     print(f"Node {node_id}:")
#     print(f"  Control Volume Area = {info['CONTROL VOLUME AREA']}")
#     print(f"  Pressure at Node = {info['PRESSURE']}")
#     for face_id, face_info in info['FACES'].items():
#         print(f"  Face {face_id}:")
#         print(f"    Midpoint = {face_info['MIDPOINT']}")
#         print(f"    Length = {face_info['LENGTH']}")
#         print(f"    Interpolated Pressure = {face_info['PRESSURE']}")
#     print()


##CREATION OF 2D AIRFOILS WITH JOUKOWSKI POISSONS FOR 2D AIRFOILS
import numpy as np
import matplotlib.pyplot as plt

def joukowsky_transform(z, c):
    return 0.2 * (z + c**2 / z)

def circle_points(radius, num_points, center):
    angles = np.linspace(0, 2 * np.pi, num_points)
    x = center[0] + radius * np.cos(angles)
    y = center[1] + radius * np.sin(angles)
    return np.column_stack((x, y))

def create_airfoil(center, radius, c, num_points):
    circle = circle_points(radius, num_points, center)
    circle_complex = circle[:, 0] + 1j * circle[:, 1]
    airfoil_complex = joukowsky_transform(circle_complex, c)
    return np.column_stack((airfoil_complex.real, airfoil_complex.imag))

def apply_transformations(airfoil, translation=(0, 0), rotation=0, scaling=(1,1)):
    # Scaling
    airfoil = airfoil * np.array(scaling)

    # Rotation
    rotation_matrix = np.array([[np.cos(rotation), -np.sin(rotation)],
                                [np.sin(rotation), np.cos(rotation)]])
    airfoil = np.dot(airfoil, rotation_matrix)

    # Translation
    airfoil = airfoil + np.array(translation)

    return airfoil

# Parameters for the rotor and stator airfoils
c_rotor = 1.0
c_stator = 1.0
circle_center_rotor = 0.1 + 0.1j
circle_center_stator = 0.1 + 0.1j
circle_radius_rotor = np.abs(c_rotor - circle_center_rotor)
circle_radius_stator = np.abs(c_stator - circle_center_stator)
num_points = 200

# Generate rotor and stator airfoils
airfoil_rotor = create_airfoil([circle_center_rotor.real, circle_center_rotor.imag], circle_radius_rotor, c_rotor, num_points)
airfoil_stator = create_airfoil([circle_center_stator.real, circle_center_stator.imag], circle_radius_stator, c_stator, num_points)

# Apply transformations to the rotor airfoil
translation_rotor = (3.4, 3.2)
rotation_rotor = np.deg2rad(-30)  # Rotate 
scaling_rotor = (2.5,-2.5)  # Scale 
airfoil_rotor_transformed = apply_transformations(airfoil_rotor, translation_rotor, rotation_rotor, scaling_rotor)

# Apply transformations to the stator airfoil
translation_stator = (1.6, 3.2)
rotation_stator = np.deg2rad(30)  # Rotate 
scaling_stator = (2.5,2.5)  # Scale 
airfoil_stator_transformed = apply_transformations(airfoil_stator, translation_stator, rotation_stator, scaling_stator)

# Plot the transformed rotor and stator airfoils
fig4 = plt.figure(figsize=(10, 5))
plt.plot(airfoil_rotor_transformed[:, 0], airfoil_rotor_transformed[:, 1], 'k-', lw=1.5, label='Rotor')
plt.plot(airfoil_stator_transformed[:, 0], airfoil_stator_transformed[:, 1], 'k-', lw=1.5, label='Stator')
plt.xlabel('LENGTH', fontsize=14)
plt.ylabel('WIDTH', fontsize=14)
#plt.title('Transformed Rotor and Stator Airfoils', fontsize=16)
plt.axis('equal')
plt.tight_layout()
#fig4.savefig('2DAIRFOIL.png', dpi=600)
plt.show()

import numpy as np
from matplotlib.path import Path
import matplotlib.colors as mcolors
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import pyamg

start = time.time()

def poisson_solver(f, dx, dy, nx, ny, num_levels):
    """
    Solve Poisson's equation on a 2D grid using the finite difference method with cell-centered scheme and multigrid.
    Assumes a zero gradient at the boundaries.
    """
    # Create the cell-centered grid
    x = np.linspace(-1+dx/2, 1-dx/2, nx)
    y = np.linspace(-1+dy/2, 1-dy/2, ny)
    X, Y = np.meshgrid(x, y, indexing='ij')
    A = pyamg.gallery.poisson((nx, ny), format='csr')

    # Create the multigrid solver
    ml = pyamg.smoothed_aggregation_solver(A, max_coarse=num_levels)
    print(ml) 

    # Solve the system
    b = dx*dy*f.flatten()
    phi = ml.solve(b, tol=1e-10, maxiter=100000)

    # Reshape the solution
    phi = phi.reshape((nx, ny))

    phi[0, :] = 0  # lower x edge
    phi[-1, :] = 0  # upper x edge
    phi[:, 0] = 0  # lower y edge
    phi[:, -1] = 0  # upper y edge

    return phi

# Define the grid size
nx, ny = 120, 120
# Define the step size of the grid
# dx, dy = 2/(nx-1), 2/(ny-1)
# # Define the grid
# x = np.linspace(-1, 1, nx)
# y = np.linspace(-1, 1, ny)
# X, Y = np.meshgrid(x, y)

# Define the step size of the grid
dx, dy = 5.25 / (nx - 1), 5 / (ny - 1)
# Define the coordinates for each dimension
x = np.linspace(-0.25, 5.25, nx)
y = np.linspace(0, 5, ny)
# Create the meshgrid
X, Y = np.meshgrid(x, y)


##MASK CREATION
# Masks for rotor and stator airfoils
def create_airfoil_mask(airfoil, X, Y):
    grid_points = np.column_stack((X.flatten(), Y.flatten()))
    airfoil_path = Path(airfoil)
    return airfoil_path.contains_points(grid_points).reshape(X.shape)

mask_rotor = create_airfoil_mask(airfoil_rotor_transformed, X, Y)
mask_stator = create_airfoil_mask(airfoil_stator_transformed, X, Y)

# Combine the masks (logical OR)
mask_combined = np.logical_or(mask_rotor, mask_stator)
#perimeter_combined = find_perimeter(mask_combined)

# Force matrix
f = np.zeros((nx, ny))

# Constant load of 100 N to the airfoil masks
force_value = -100
f[mask_combined] = force_value
#total_force = np.sum(f[mask_combined])

# Boundary points of the airfoils
boundary_rotor = np.vstack((airfoil_rotor_transformed[:-1], airfoil_rotor_transformed[1:]))
boundary_stator = np.vstack((airfoil_stator_transformed[:-1], airfoil_stator_transformed[1:]))

# New array to store the opposite forces
boundary_mask = np.zeros_like(X)

# Application of opposite force along the boundary lines and update the boundary_mask array
def apply_opposite_force_boundary(boundary, X, Y, f, boundary_mask, force_value):
    for i in range(len(boundary) - 1):
        point_a, point_b = boundary[i], boundary[i + 1]
        idx_a = np.argmin(np.linalg.norm(np.column_stack((X.flatten() - point_a[0], Y.flatten() - point_a[1])), axis=1))
        idx_b = np.argmin(np.linalg.norm(np.column_stack((X.flatten() - point_b[0], Y.flatten() - point_b[1])), axis=1))
        x_idx_a, y_idx_a = np.unravel_index(idx_a, X.shape)
        x_idx_b, y_idx_b = np.unravel_index(idx_b, X.shape)
        f[x_idx_a, y_idx_a] += -force_value
        f[x_idx_b, y_idx_b] += -force_value
        boundary_mask[x_idx_a, y_idx_a] = -force_value
        boundary_mask[x_idx_b, y_idx_b] = -force_value

# Apply the opposite force along the boundary lines of the airfoils
apply_opposite_force_boundary(boundary_rotor, X, Y, f, boundary_mask, -force_value)
apply_opposite_force_boundary(boundary_stator, X, Y, f, boundary_mask, -force_value)

# force_perimeter = 50
# f[perimeter_combined] = force_perimeter

f[0,:] = 0
f[-1,:] = 0
f[:,0] = 0
f[:,-1] = 0

# Solve Poisson's equation with multigrid
num_levels = 16
phi = poisson_solver(f, dx, dy, nx, ny, num_levels)

print(f'duration: {time.time() - start}')

# 3D plot of the solution
fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection='3d')
ax2.plot_surface(X, Y, phi)
#ax2.set_title("3D Plot")
#fig2.savefig("3D_PHIPlot.svg", format='svg', dpi=1200)
plt.show()

# 2D plot of the solution
fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
im = ax3.imshow(phi, origin='lower', extent=(-1, 1, -1, 1), cmap='viridis')
plt.colorbar(im, ax=ax3)
#ax3.set_title("2D Plot")
#ax3.set_xlabel('LENTGH')
#ax3.set_ylabel('WIDTH')
#fig3.savefig("2D_PHIPlot.png", dpi=600)
#fig3.savefig("2D_contour.svg", format='svg', dpi=1200)

plt.show()


##INTERPOLATION OF PHI IN 2D GRID PATTERN 
phi_interpolator = interpolate.interp2d(X, Y, phi, kind='linear')

# Phi values at the midpoints of the faces
for node_id, info in node_info.items():
    for face_id, face_info in info['FACES'].items():
        # Calculate the interpolated phi value
        interpolated_phi = phi_interpolator(*face_info['MIDPOINT'])

        # Store the interpolated phi value
        face_info['PHI'] = interpolated_phi

# # Print the interpolated phi values for each node
# for node_id, info in node_info.items():
#     print(f"Node {node_id}:")
#     for face_id, face_info in info['FACES'].items():
#         print(f"  Face {face_id}: Interpolated Phi Value = {face_info['PHI']}")


# Calculate the final pressure values at the midpoints of the faces
for node_id, info in node_info.items():
    for face_id, face_info in info['FACES'].items():
        # Calculate the final pressure value
        final_pressure = face_info['PHI'] + face_info['PRESSURE']
        # Store the final pressure value
        face_info['FINAL PRESSURE'] = final_pressure

# # DETAILS FOR EACH NODE
# for node_id, info in node_info.items():
#     print(f"Node {node_id}:")
#     print(f"  Control Volume Area = {info['CONTROL VOLUME AREA']}")
#     print(f"  Pressure at Node = {info['PRESSURE']}")
#     for face_id, face_info in info['FACES'].items():
#         print(f"  Face {face_id}:")
#         print(f"    Midpoint = {face_info['MIDPOINT']}")
#         print(f"    Length = {face_info['LENGTH']}")
#         print(f"    Interpolated Phi = {face_info['PHI']}")
#         print(f"    Interpolated Pressure = {face_info['PRESSURE']}")
#         print(f"    Final Pressure = {face_info['FINAL PRESSURE']}")
#     print()


# Update the face details
for node_id, info in node_info.items():
    for face_id, face_info in info['FACES'].items():
        centroid = np.array(list(face_id)[0])
        midpoint = np.array(face_info['MIDPOINT'])

        # Calculate the normal vector
        normal = midpoint - centroid
        normal /= np.linalg.norm(normal)  # Normalize the normal vector

        face_info['NORMAL'] = normal

# Calculate the forces on each node
for node_id, info in node_info.items():
    total_force = np.zeros(2)
    for face_id, face_info in info['FACES'].items():
        force = np.array(face_info['FINAL PRESSURE']) * np.array(face_info['NORMAL']) * face_info['LENGTH']
        total_force += force
    info['TOTAL FORCE'] = total_force  # Add the total force to the node info

# Velocities for each node
for node_id, info in node_info.items():
    info['VELOCITY'] = np.zeros(2)

# Update the velocities
dt = 0.025  # Time step
for node_id, info in node_info.items():
    acceleration = info['TOTAL FORCE'] / mass
    info['VELOCITY'] += acceleration * dt


# Print the details for each node
for node_id, info in node_info.items():
    print(f"Node {node_id}:")
    print(f"  Control Volume Area = {info['CONTROL VOLUME AREA']}")
    print(f"  Pressure at Node = {info['PRESSURE']}")
    print(f"    total_force = {info['TOTAL FORCE']}")
    print(f"  Velocity = {info['VELOCITY']}")

    for face_id, face_info in info['FACES'].items():
        print(f"  Face {face_id}:")
        print(f"    Midpoint = {face_info['MIDPOINT']}")
        print(f"    Length = {face_info['LENGTH']}")
        print(f"    Interpolated Phi = {face_info['PHI']}")
        print(f"    Interpolated Pressure = {face_info['PRESSURE']}")
        print(f"    Final Pressure = {face_info['FINAL PRESSURE']}")
    print()


# Update the positions of the nodes
for node_id, info in node_info.items():
    node_coords[node_id] += info['VELOCITY'] * dt
# Calculate the time step for each node
for node_id, info in node_info.items():
    velocity_magnitude = np.linalg.norm(info['VELOCITY'])
    info['TIME STEP'] = info['CONTROL VOLUME AREA'] / velocity_magnitude
# Calculate the global time step
time_steps = [info['TIME STEP'] for info in node_info.values()]
global_time_step = min(time_steps)

# Print the updated positions for each node
for node_id, info in node_info.items():
    print(f"Node {node_id}: New Position = {node_coords[node_id]}")

# Plot the new positions of the nodes
fig1 = plt.figure(figsize=(10, 10))
plt.scatter(node_coords[:, 0], node_coords[:, 1])

# Draw connections between the nodes
for connection in pattern_connections:
    node1 = node_coords[connection[0]]
    node2 = node_coords[connection[1]]
    plt.plot([node1[0], node2[0]], [node1[1], node2[1]], 'k-')

plt.xlabel("x")
plt.ylabel("y")
plt.axis('equal')
plt.grid(True)
fig1.savefig("final.svg", format='svg', dpi=1200)
plt.show()

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
