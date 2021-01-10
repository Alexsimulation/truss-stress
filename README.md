### TRUSS STRESS


# What it is
A 3D truss stress analysis tool. Returns the stress in each edge of the truss, given geometry and loads. 


# How to run
Two versions of the program are provided:
- **TS** is a user-oriented version that opens a command window when executed, asks for the truss file name, prints the stresses in the command window and in a file named *stress.csv*. The user can't choose another output file. The user can also execute a truss file with the TS application.
- **TSB** is a command line executed version that automatically exits execution when it's done. To run it, open a command window in the directory of the TSB application, and write the command *TSB filename.frm outfilename.csv*, the output file being optional. This version allows user to choose a custom name for the output file.


# Truss definition file
Use a custom made filename.trs text file to define the truss. To write the truss file, use the following tags:

- **// START TRUSS //** : Start of the truss text file.
- **E: value** : Defines the young's modulus of the material, in Pa.
- **d: value** : Defines the density of the material, in kg/m^3.
- **v: value_x, value_y, value_z** : Defines a vertex with its position, in meters.
- **f: v_index ; value_x, value_y, value_z** : Defines an external force on a vertex, in newtons.
- **c: v_index** : Constrains a predefined vertex. Constrained vertex are fixed in position.
- **e: area, v0_index, v1_index** : Defines an edge with its cross-sectional area, its starting vertex and its end vertex.
- **// END TRUSS //** : End of the truss text file.

Vertex and edges are indexed in line order, starting at index 0. An example truss file is provided, named *example.trs*.


# Version and credits
Version : 0.1 - private release
(c) 2021 Alexis Angers (https://github.com/Alexsimulation). Private and educational use only.