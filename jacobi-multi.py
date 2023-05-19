##JACOBI AN DMULTIGRID METHODS COMPARISONS WITH ERROR DISTRIBUTION AND CONVERGENCE MAP

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
from functools import partial

start = time.time()

def poisson_solver(f, dx, dy, nx, ny):
    """
    Solve Poisson's equation on a 2D grid using the Jacobi method.
    Assumes a zero gradient at the boundaries.
    """
    # Initialize the solution array with zeros
    phi = np.zeros((nx, ny))
    # Set the initial condition
    #phi[0,0] = 0.
    # Set the boundary conditions
    phi[-1, :] = 0.
    phi[:, -1] = 0.
    phi[0, :] = 0.
    phi[:, 0] = 0.
    # Initialize the L2 norm error
    error = 1.0
    error_historyj = []
    # Iterate until the solution converges or until the maximum number of iterations is reached
    max_iterations = 2000
    for i in range(max_iterations):
        # Update the solution at every point using the Jacobi method
        phi_new = (phi[:-2,1:-1] + phi[2:,1:-1] + phi[1:-1,:-2] + phi[1:-1,2:])/4. - f[1:-1,1:-1]*dx*dy/4.
        # Compute the L2 norm error
        error = np.sqrt(np.sum((phi_new - phi[1:-1,1:-1])**2))
        # Update the solution
        phi = np.pad(phi_new, ((1,1), (1,1)), 'constant')
        # Stop the iteration if the L2 norm error falls below the threshold
        if error < 1e-8:
            break
        error_historyj.append(error)
    return phi, error_historyj

# Define the grid size
nx, ny = 100, 100
# Define the step size of the grid
dx, dy = 2/(nx-1), 2/(ny-1)
# Define the grid
x = np.linspace(-1, 1, nx)
y = np.linspace(-1, 1, ny)
X, Y = np.meshgrid(x, y)
# Define the function f and change magnitude 
f = 1000*np.ones((nx, ny))
f[0,:] = 0
f[-1,:] = 0
f[:,0] = 0
f[:,-1] = 0
f = f*100/((nx-2)*(ny-2))

# Solve Poisson's equation
phi_jacobi, error_history_jacobi = poisson_solver(f, dx, dy, nx, ny)

# Visualize the solution
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, phi_jacobi)

print(f'duration: {time.time() - start}')

print(phi_jacobi)

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

    # Solve the system
    b = dx*dy*f.flatten()
    residual_history = []  # Initialize an empty list for the residuals
    phi = ml.solve(b, tol=1e-10, maxiter=20000, accel='cg', residuals=residual_history)


    # Reshape the solution
    phi = phi.reshape((nx, ny))

    phi[0, :] = 0  # lower x edge
    phi[-1, :] = 0  # upper x edge
    # Neumann boundary conditions
    # phi[0, :] = phi[1, :]  # lower x edge
    # phi[-1, :] = phi[-2, :]  # upper x edge
    phi[:, 0] = 0  # lower y edge
    phi[:, -1] = 0  # upper y edge

    return phi, residual_history

# Define the grid size
nx, ny = 100, 100
# Define the step size of the grid
dx, dy = 2/(nx-1), 2/(ny-1)
# Define the grid
x = np.linspace(-1, 1, nx)
y = np.linspace(-1, 1, ny)
X, Y = np.meshgrid(x, y)
# Define the function f and change magnitude 
f = -1000*np.ones((nx, ny))
f[0,:] = 0
f[-1,:] = 0
f[:,0] = 0
f[:,-1] = 0
f = f*100/((nx-2)*(ny-2))

# Solve Poisson's equation with multigrid
num_levels = 16
phi_multigrid, residual_history_multigrid = poisson_solver(f, dx, dy, nx, ny, num_levels)

# Visualize the solution
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, phi_multigrid)
#fig.savefig(".svg", format='svg', dpi=1200)

# fig1 = plt.figure()
# ax1 = fig1.add_subplot(111, projection='3d')
# ax1.plot_surface(X_def, Y_def, phi)

print(f'duration: {time.time() - start}')

print(phi_multigrid)

fig = plt.figure()
plt.plot(error_history_jacobi, label='Jacobi')
plt.plot(residual_history_multigrid, label='Multigrid')
plt.xlabel('Iteration')
plt.ylabel('Convergence history of error')
plt.yscale('log')
plt.legend()
plt.grid()
#fig.savefig("errorhistory.svg", format='svg', dpi=1200)

# Compute the error between the two solutions
error = abs(phi_jacobi - phi_multigrid)

# Create a heatmap of the error distribution
fig = plt.figure()
plt.imshow(error, extent=(-1, 1, -1, 1), origin='lower', cmap='viridis')
plt.colorbar(label='Error')
plt.xlabel('x')
plt.ylabel('y')
fig.savefig("error.svg", format='svg', dpi=1200)

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

