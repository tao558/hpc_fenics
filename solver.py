from __future__ import print_function # import the print functions
from fenics import * # import fenics functions
from dolfin import * # import fenics functions
from mshr import * # import mshr functions
from math import sin, cos, pi # import sin, cos, and pi
import matplotlib.pyplot as plt

n = 10    # number of windings
mesh = Mesh('mesh.xml')



# Define function space from the mesh
V = FunctionSpace(mesh, 'P', 1)


# Define boundary condition
bc = DirichletBC(V, Constant(0), 'on_boundary')


markers = MeshFunction('size_t', mesh, 2, mesh.domains()) # Define subdomain markers
dx = Measure('dx', domain=mesh, subdomain_data=markers) # Redefine integration measure


# Define northern current densities +1 A
J_N = Constant(1.0)
# Define southern current densities -1 A
J_S = Constant(-1.0)


# Define magnetic permeability
class Permeability(UserExpression):
# constructor for permeability class
    def __init__(self, markers, **kwargs):
        super().__init__(**kwargs)
        self.markers = markers
# define function to evaluate each cell
    def eval_cell(self, values, x, cell):
#if the cell is flagged zero, it’s a vacuum, so use this value
        if self.markers[cell.index] == 0:
            values[0] = 4*pi*1e-7 # vacuum
#if the cell is flagged one, it’s iron, so use this value
        elif self.markers[cell.index] == 1:
            values[0] = 1e-5      # iron (should really be 6.3e-3)
# otherwise the cell is copper
        else:
            values[0] = 1.26e-6   # copper


mu = Permeability(markers, degree=1) # Define mu using iron (degree=1)


# Define variational problem
A_z = TrialFunction(V) # Create a trialfunction argument
v = TestFunction(V) # Create a testfunction argument
a = (1 / mu)*dot(grad(A_z), grad(v))*dx # Compute LHS of equation
L_N = sum(J_N*v*dx(i) for i in range(2, 2 + n)) # Set up northward current density integral over mesh
L_S = sum(J_S*v*dx(i) for i in range(2 + n, 2 + 2*n)) # Set up southward current density integral over mesh
L = L_N + L_S # Sum of northward and southward


# Solve variational problem
A_z = Function(V) # Declare variable to store potential
solve(a == L, A_z, bc) # Calculate potential by solving equation and store in A_z. Apply BC


# Compute magnetic field (B = curl A)
W = VectorFunctionSpace(mesh, 'P', 1) # Result vector
B = project(as_vector((A_z.dx(1), -A_z.dx(0))), W) # Compute curl

if MPI.rank(MPI.comm_world) == 0:
    print("Got here")
    # Plot solution
    plot(A_z) # Plot of the z-component Az of the magnetic vector potential.
    plot(B) # Plot of the magnetic field B in the xy-plane.

    # Hold plot
    # plt.show()

# Save solution to file
vtkfile_A_z = File('magnetostatics/potential.pvd') # Create file for potential
vtkfile_B = File('magnetostatics/field.pvd') # Create file for magnetic field

vtkfile_A_z << A_z # Write out potential
vtkfile_B << B # Write out magnetic field

