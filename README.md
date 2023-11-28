# Hydrosolver2D
This code is meant to solve the 2D euler equations for arbitrary starting distributions and periodic boundary conditions. It uses a dimensionally split finite volume approach on a fixed grid and solves advection and source terms seperatly. 

For the advection step one can choose different flux functions implemented in flux.h. The equation used for determining pressure in the source terms can be chosen from pressure.h.
