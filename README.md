# 4D-Mesh
Third Year Individual Research Project (IRP) concerning the realisation of a Computational Fluid Dynamics (CFD) algorithm for the simulation of a rotor and stator in a 4D mesh.

Advanced 4D Mesh Generation for Rotor-Stator Airfoils in Computational Fluid Dynamics

The project essentially consists of creating a body fitting mesh algorithm in python to simulate the interaction between a rotor and stator moving at relatively high speeds. The algorithm was built for a 2D, 3D and 4D model. It uses as a solid basis an Arbitrary Lagrangian-Eulerian (ALE) method which involves the creation of the rotor and stator airfoils, the implementation of Poisson's equation for grid deformation in 2D, 3D, and 4D cases, the use of a new mesh generation technique and the derivation of control volumes, interface pressures and updated nodes position.
The project still has several limitations due to the quality of the code and the way the volumes of the 3-dimensional faces of the control volumes in the 4D mesh are calculated. Subsequently, if this problem will be solved, the simulation of the model over time will take over. In this repository the project and the essential codes will be found. 
