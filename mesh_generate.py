from __future__ import print_function # import the print functions
from fenics import * # import fenics functions
from mshr import * # import mshr functions
from math import sin, cos, pi # import sin, cos, and pi
import matplotlib.pyplot as plt


a = 1.0   # inner radius of iron cylinder
b = 1.2   # outer radius of iron cylinder
c_1 = 0.8 # radius for inner circle of copper wires
c_2 = 1.4 # radius for outer circle of copper wires
r = 0.1   # radius of copper wires
R = 5.0   # radius of domain
n = 10    # number of windings


domain = Circle(Point(0, 0), R) # Define geometry for background


cylinder = Circle(Point(0, 0), b) - Circle(Point(0, 0), a) # Define geometry for iron cylinder


# Define geometry for wires (N = North (up), S = South (down))
angles_N = [i*2*pi/n for i in range(n)] # define an array of northern angles
angles_S = [(i + 0.5)*2*pi/n for i in range(n)] # define an array of southern angles


# define circle of wires for each northern angle
wires_N = [Circle(Point(c_1*cos(v), c_1*sin(v)), r) for v in angles_N]


# define circle of wires for each southern angle
wires_S = [Circle(Point(c_2*cos(v), c_2*sin(v)), r) for v in angles_S]


domain.set_subdomain(1, cylinder) # Set subdomain for iron cylinder


# Set subdomains for all northern wires
for (i, wire) in enumerate(wires_N):
    domain.set_subdomain(2 + i, wire)


# Set subdomain for all southern wires
for (i, wire) in enumerate(wires_S):
    domain.set_subdomain(2 + n + i, wire)


# Create mesh with granularity of 128
mesh = generate_mesh(domain, 128)

File("mesh.xml") << mesh