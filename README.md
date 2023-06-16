# BlochQT 
Simulation program that couples Bloch wave image simulations with quantum
trajectory calculations to provide a hydrodynamic interpretation to images
obtained in transmission electron microscopy. The details of the simulation are
provided in the following publication: 

<p>Rudinsky, S & Gauvin, R.<br>
A hydrodynamic approach to electron beam imaging using a Bloch wave representation
Ultramicroscopy, 2020, 212, 112979<br>
doi: (https://doi.org/10.1016/j.ultramic.2020.112979)</p>

# Build
Project requires Boost, Eigen, and FFTW3
## Linux
Builds with cmake. It is recommended to do the following:

```
mkdir build
cd build/
cmake ..
make
```

Binaries will be created in bin/ directory from project

## Windows
Project can also be built in VS using CMakeLists.txt file 

# Simulation inputs
## Material input file
For the material file, program accepts a .cfg
Here is an example, lines that begin with a hash are comments. 

```
Crystal Type = FCC, BCC, DC, or NONE 
# First vector of unit call
H0(1,1) = 3.6149 
H0(1,2) = 0
H0(1,3) = 0
H0(2,1) = 0
# Second vector of unit cell
H0(2,2) = 3.6149 
H0(2,3) = 0
# Third vector of unit cell
H0(3,1) = 0
H0(3,2) = 0
H0(3,3) = 3.6149 
Number Of Elements = 2

# Elements in model and associated Debye-Waller factors

Cu
17  
Ni
3

# Number of atoms in unit call, coordinates of each atom and which element type it is
# The element type is indicated by the order with which they are specified in the
# previous section.

Number Of Atoms = 4
0, 0, 0, 1 
0.5, 0.5, 0, 2
0, 0.5, 0.5, 1
0.5, 0, 0.5, 1
```

## Program inputs
To create an input file for Planewave the order is the following:

- Energy in eV
- Name of material file including location
- Tilt? (y/n)
- If y then choose Bragg (B), angle (A) or systematic row (S)
- Input Bragg reflection or angle depending on choice
- If systematic row, then input beam and integer tilt
- Number of grid points in x and y (xg yg)
- Absorption? (y/n)
- Strong and weak beam limit (c1 c2)
- Max thickness in angstrom
- Number of slices
- Compute trajectories (y/n)
- Number of trajectories (nx ny)

# TODO
- Probe simulations produce incorrect results.
