import plotly.graph_objs as go
import pyvista as pv
import numpy as np
import math
from numpy.linalg import det
import itertools
import matplotlib.pyplot as plt

# Leggi i dati dal file .txt
##CHANGE IMPORTED FILE NAME IF YOU WANT TO SEE ORIGINAL DISCONNECTED NODES.
with open('MeshPatch4D_SSK3.txt', 'r') as f:
    lines = f.readlines()
# Salva i dati in due liste separate
simplices = [list(map(int, l.strip().split())) for l in lines[5:3701]]
nodes = [list(map(float, l.strip().split())) for l in lines[3705:4176]]
# Crea un dizionario per il mesh 4D
mesh = {}
# Crea una lista di nodi vuota
nodes_list = [None] * len(nodes)
# Aggiungi i nodi al dizionario e alla lista di nodi
for node in nodes:
    node_id = int(node[0])
    coords = tuple(node[1:])
    mesh[node_id] = {'coords': coords}
    nodes_list[node_id - 1] = node  #lista di tutti i nodi con indici file

# Aggiungi i simplessi al dizionario
for simplex in simplices:
    simplex_id = simplex[0]
    node_ids = simplex[1:]
    tet_coords = tuple(mesh[node_id]['coords'] for node_id in node_ids)
    mesh[simplex_id] = {'type': 'tetrahedron', 'nodes': node_ids, 'coords': tet_coords}

# Esempio di come accedere alle coordinate di un nodo e le coordinate dei nodi di un simplex
simplex_id = 3553
simplex_node_ids = mesh[simplex_id]['nodes']
simplex_node_coords = [nodes_list[node_id - 1] for node_id in simplex_node_ids]
print(f"I nodi del simplex con id {simplex_id} sono: {simplex_node_ids}")
print(f"Le coordinate dei nodi del simplex con id {simplex_id} sono:")
for node in simplex_node_coords:
    print(node)

def find_disconnected_nodes(simplices, num_nodes):
    connected_nodes = set()
    for simplex in simplices:
        for node in simplex[1:]:
            connected_nodes.add(node)
    disconnected_nodes = set(range(1, num_nodes + 1)) - connected_nodes
    return disconnected_nodes

##disocnnected nodes
# 82, 83, 86, 88, 98, 99, 106, 107, 109, 111, 124, 125, 136, 137, 139, 140, 143, 144, 153, 154, 155, 157, 161, 162

import plotly.graph_objs as go
import numpy as np
from plotly.subplots import make_subplots

# def visualize_4d_mesh(nodes_list, simplices, disconnected_nodes, chosen_node=470):
#     unique_4th_coords = sorted(set(node[4] for node in nodes_list))

#     num_subplots = len(unique_4th_coords)
#     subplot_titles = [f"4th Coord: {coord}" for coord in unique_4th_coords]

#     fig = make_subplots(rows=1, cols=num_subplots, specs=[[{'type': 'scatter3d'}] * num_subplots], subplot_titles=subplot_titles)

#     for idx, fourth_coord in enumerate(unique_4th_coords):
#         # Filter nodes with the current 4th coordinate value
#         filtered_nodes = [node for node in nodes_list if node[4] == fourth_coord]
#         filtered_coords = np.array([node[1:4] for node in filtered_nodes])

#         # Set colors for disconnected nodes
#         node_colors = ['red' if node[0] in disconnected_nodes else 'blue' for node in filtered_nodes]

#         # Plot the nodes
#         fig.add_trace(go.Scatter3d(x=filtered_coords[:, 0], y=filtered_coords[:, 1], z=filtered_coords[:, 2],
#                                    mode='markers',
#                                    marker=dict(color=node_colors, size=5),
#                                    text=[f'Node ID: {node[0]}<br>Coords: {node[1:4]}<br>4th Coord: {node[4]}' for node in filtered_nodes],
#                                    hoverinfo='text',
#                                    name=f'4th Coord: {fourth_coord}',
#                                    legendgroup=f'group{idx}',
#                                    showlegend=False),
#                      row=1, col=idx+1)

#         # Draw connections (simplices) around the chosen node
#         if chosen_node:
#             for simplex in simplices:
#                 if chosen_node in simplex[1:]:
#                     lines = itertools.combinations(simplex[1:], 2)
#                     for line in lines:
#                         node1 = nodes_list[line[0] - 1]
#                         node2 = nodes_list[line[1] - 1]
#                         if node1[4] == node2[4] == fourth_coord:
#                             fig.add_trace(go.Scatter3d(x=[node1[1], node2[1]], y=[node1[2], node2[2]], z=[node1[3], node2[3]],
#                                                        mode='lines',
#                                                        line=dict(color='green', width=2),
#                                                        hoverinfo='none',
#                                                        showlegend=False),
#                                          row=1, col=idx+1)

#     # Add labels to the plot
#     fig.update_layout(scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Z'),
#                      title='Separate 3D Spaces for Each Unique 4th Coordinate Value',width=1200, height=1200)

#     # Show the plot
#     fig.show()

def visualise_4d_mesh_separate_windows(nodes_list, simplices, disconnected_nodes, chosen_node=325):
    unique_4th_coords = sorted(set(node[4] for node in nodes_list))

    for idx, fourth_coord in enumerate(unique_4th_coords):
        fig = go.Figure()

        # Filter nodes with the current 4th coordinate value
        filtered_nodes = [node for node in nodes_list if node[4] == fourth_coord]
        filtered_coords = np.array([node[1:4] for node in filtered_nodes])

        # Set colors for disconnected nodes
        node_colors = ['red' if node[0] in disconnected_nodes else 'blue' for node in filtered_nodes]

        # Plot the nodes
        fig.add_trace(go.Scatter3d(x=filtered_coords[:, 0], y=filtered_coords[:, 1], z=filtered_coords[:, 2],
                                   mode='markers',
                                   marker=dict(color=node_colors, size=5),
                                   text=[f'Node ID: {node[0]}<br>Coords: {node[1:4]}<br>4th Coord: {node[4]}' for node in filtered_nodes],
                                   hoverinfo='text',
                                   name=f'4th Coord: {fourth_coord}'))

        # Draw connections (simplices) around the chosen node
        if chosen_node:
            for simplex in simplices:
                if chosen_node in simplex[1:]:
                    lines = itertools.combinations(simplex[1:], 2)
                    for line in lines:
                        node1 = nodes_list[line[0] - 1]
                        node2 = nodes_list[line[1] - 1]
                        if node1[4] == node2[4] == fourth_coord:
                            fig.add_trace(go.Scatter3d(x=[node1[1], node2[1]], y=[node1[2], node2[2]], z=[node1[3], node2[3]],
                                                       mode='lines',
                                                       line=dict(color='green', width=2),
                                                       hoverinfo='none',
                                                       showlegend=False))

        # Add labels to the plot
        fig.update_layout(scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Z'),
                         title=f'Separate 3D Space for 4th Coordinate Value: {fourth_coord}',width=1200, height=800)

        # Show the plot
        fig.show()


# Find disconnected nodes
num_nodes = len(nodes)
disconnected_nodes = find_disconnected_nodes(simplices, num_nodes)

# # Visualize the 4D mesh with disconnected nodes highlighted
# visualize_4d_mesh(nodes_list, simplices, disconnected_nodes)
#visualize_4d_mesh(nodes_list, simplices, disconnected_nodes, chosen_node=470)
visualise_4d_mesh_separate_windows(nodes_list, simplices, disconnected_nodes, chosen_node=325)


def find_simplices_for_node(chosen_node, simplices, nodes_list):
    simplices_in_same_4d = []
    simplices_between_4d = []
    for simplex in simplices:
        if chosen_node in simplex[1:]:
            # Check if all the nodes in the simplex have the same 4th coordinate
            fourth_coords = [nodes_list[node_id - 1][4] for node_id in simplex[1:]]
            if len(set(fourth_coords)) == 1:
                simplices_in_same_4d.append(simplex)
            else:
                simplices_between_4d.append(simplex)

    return simplices_in_same_4d, simplices_between_4d
chosen_node = 470
simplices_in_same_4d, simplices_between_4d = find_simplices_for_node(chosen_node, simplices, nodes_list)

print("Simplices in the same 4D system:")
for simplex in simplices_in_same_4d:
    print(simplex)

print("\nSimplices between different 4D systems:")
for simplex in simplices_between_4d:
    print(simplex)



# def check_all_simplices(simplices, nodes_list):
#     for simplex in simplices:
#         # Get the 4th coordinate for each node in the simplex
#         fourth_coords = [nodes_list[node_id - 1][4] for node_id in simplex[1:]]
        
#         # If all the 4th coordinates are the same, then return False
#         if len(set(fourth_coords)) == 1:
#             return False
#     return True
# all_simplices_have_different_4d = check_all_simplices(simplices, nodes_list)
# print(f"All simplices have nodes connected between different 4D systems: {all_simplices_have_different_4d}")



# from collections import defaultdict

# def find_overlapping_nodes(nodes):
#     coord_groups = defaultdict(list)
#     for node in nodes:
#         node_id = int(node[0])
#         coord_3d = tuple(node[1:4])
#         coord_groups[coord_3d].append(node_id)

#     overlapping_nodes = {}
#     for coords, node_ids in coord_groups.items():
#         if len(node_ids) > 1:
#             overlapping_nodes[coords] = node_ids

#     return overlapping_nodes

# overlapping_nodes = find_overlapping_nodes(nodes)
# print("Overlapping nodes:")
# for coords, node_ids in overlapping_nodes.items():
#     print(f"Coords {coords}: Node IDs {node_ids}")













































# import itertools

# def visualize_single_4d_mesh(nodes_list, simplices, disconnected_nodes, chosen_fourth_coord):
#     filtered_nodes = [node for node in nodes_list if node[4] == chosen_fourth_coord]
#     filtered_coords = np.array([node[1:4] for node in filtered_nodes])

#     node_colors = ['red' if node[0] in disconnected_nodes else 'blue' for node in filtered_nodes]

#     fig = go.Figure()

#     fig.add_trace(go.Scatter3d(x=filtered_coords[:, 0], y=filtered_coords[:, 1], z=filtered_coords[:, 2],
#                                mode='markers',
#                                marker=dict(color=node_colors, size=5),
#                                text=[f'Node ID: {node[0]}<br>Coords: {node[1:4]}<br>4th Coord: {node[4]}' for node in filtered_nodes],
#                                hoverinfo='text',
#                                name=f'4th Coord: {chosen_fourth_coord}'))

#     for simplex in simplices:
#         lines = itertools.combinations(simplex[1:], 2)
#         for line in lines:
#             node1 = nodes_list[line[0] - 1]
#             node2 = nodes_list[line[1] - 1]
#             if node1[4] == node2[4] == chosen_fourth_coord:
#                 fig.add_trace(go.Scatter3d(x=[node1[1], node2[1]], y=[node1[2], node2[2]], z=[node1[3], node2[3]],
#                                            mode='lines',
#                                            line=dict(color='green', width=2),
#                                            hoverinfo='none',
#                                            showlegend=False))

#     fig.update_layout(scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Z'),
#                      title=f'3D Mesh for 4th Coordinate Value {chosen_fourth_coord}')

#     fig.show()

# # Visualize the first 4D system (4th coordinate = 0)
# visualize_single_4d_mesh(nodes_list, simplices, disconnected_nodes, chosen_fourth_coord=0)



# import itertools

# def visualize_4d_mesh(nodes_list, simplices, disconnected_nodes):
#     unique_4th_coords = sorted(set(node[4] for node in nodes_list))

#     num_subplots = len(unique_4th_coords)
#     subplot_titles = [f"4th Coord: {coord}" for coord in unique_4th_coords]

#     fig = make_subplots(rows=1, cols=num_subplots, specs=[[{'type': 'scatter3d'}] * num_subplots], subplot_titles=subplot_titles)

#     for idx, fourth_coord in enumerate(unique_4th_coords):
#         # Filter nodes with the current 4th coordinate value
#         filtered_nodes = [node for node in nodes_list if node[4] == fourth_coord]
#         filtered_coords = np.array([node[1:4] for node in filtered_nodes])

#         # Set colors for disconnected nodes
#         node_colors = ['red' if node[0] in disconnected_nodes else 'blue' for node in filtered_nodes]

#         # Plot the nodes
#         fig.add_trace(go.Scatter3d(x=filtered_coords[:, 0], y=filtered_coords[:, 1], z=filtered_coords[:, 2],
#                                    mode='markers',
#                                    marker=dict(color=node_colors, size=5),
#                                    text=[f'Node ID: {node[0]}<br>Coords: {node[1:4]}<br>4th Coord: {node[4]}' for node in filtered_nodes],
#                                    hoverinfo='text',
#                                    name=f'4th Coord: {fourth_coord}',
#                                    legendgroup=f'group{idx}',
#                                    showlegend=False),
#                      row=1, col=idx+1)

#         # Draw tetrahedrons for the current 4D system
#         for simplex in simplices:
#             lines = itertools.combinations(simplex[1:], 2)
#             for line in lines:
#                 node1 = nodes_list[line[0] - 1]
#                 node2 = nodes_list[line[1] - 1]
#                 if node1[4] == node2[4] == fourth_coord:
#                     fig.add_trace(go.Scatter3d(x=[node1[1], node2[1]], y=[node1[2], node2[2]], z=[node1[3], node2[3]],
#                                                mode='lines',
#                                                line=dict(color='green', width=2),
#                                                hoverinfo='none',
#                                                showlegend=False),
#                                  row=1, col=idx+1)

#     # Add labels to the plot
#     fig.update_layout(scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Z'),
#                      title='Separate 3D Spaces for Each Unique 4th Coordinate Value')

#     # Show the plot
#     fig.show()

# # Find disconnected nodes
# num_nodes = len(nodes)
# disconnected_nodes = find_disconnected_nodes(simplices, num_nodes)

# # Visualize the 4D mesh with disconnected nodes highlighted and tetrahedrons drawn
# visualize_4d_mesh(nodes_list, simplices, disconnected_nodes)



