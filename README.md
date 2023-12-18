# Hydrosolver2D
This code solves the 2D euler equations for arbitrary starting distributions and periodic boundary conditions. It uses a dimensionally split finite volume approach on a fixed grid and solves advection and source terms seperatly. 

For the advection step one can choose different flux functions implemented in flux.h. The equation used for determining pressure in the source terms can be chosen from pressure.h.

The simulation itself can be run manually from the main.cpp file, MakeFile support is still to be added. It results in a txt file which then can be visualized via the plotGrid.py script.
