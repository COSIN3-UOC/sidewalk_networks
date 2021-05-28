# OD-adjusted weighted edge betweenness for sidewalk networks

The code in this repository can be used to calculate the edge (and node) betweenness of sidewalk networks using empirically-grounded OD matrices estimated from geographical population and point-of-interest (POI) data. The function to estimate the OD matrix and to calculate the network betweenness is a user-defined MATLAB MEX function implemented in C++. Please cite [1] in any works using or derived from this software.
 
## System Requirements 	

* Matlab 2019b or newer
* C++ Boost library version 1.64.0 or newer (https://www.boost.org/)   
        
## Inputs:

The main function, ```get_ebw_make_od```, takes 6 inputs:
1) a MATLAB Sparse Matrix with edge lengths (in meters, feet, etc.) as weights
2) N, the number of nodes in the network
3) the number of layers in the network. This should be set to 1, as multilayer functionality still needs to be implemented.
4) a 1 x N array of integers indicating the number of POIs assigned to each network node
4) a 1 x N array of doubles indicating the population (i.e. number of pedestrians) assigned to each network node
6) a 2 x N array of doubles indicating the x and y coordinates (respectively) of each network node

and returns 2 outputs:
1) a 1 x N array of doubles indicating the node betweenness of each nod
2) an N x N array of doubles indicating the edge betweenness of each directed edge (i,j)

An example call from within MATLAB looks like:
```
[NB, directed_EB] = get_ebw_make_od(A_sparse, N, 1, vStores, vPop, vCoords);
```
An example MATLAB script with data is also provided in this repository, as described below.

## Executing the code

First, download the contents of this repository. 

### Compiling the MEX function

The MEX function must be compiled before use. In the MATLAB command line, move to the directory where the repository is stored 

```
cd /path/to/directory/ ;
```

Next, compile the function using the ```mex``` command. Make sure to write the appropriate path to the C++ Boost library on your machine. If you have a newer version than 1.64.0, change that part of the path below to the correct version.

```
mex -largeArrayDims get_ebw_make_od.cpp -I/path/to/boost/1.64.0/include/ ;
```

NOTE: To improve the performance of the software, some data structures used are defined at compile time and thus remain static. Accordingly, if networks larger than 65,000 nodes are used, the constant MAX_NODES defined on line 16 should be changed to the appropriate number. Likewise, MAX_COLA (line 15) should be changed to MAX_NODES + 1. The script will then need to be recompiled. 

### Using the example script

An example MATLAB script (```get_directed_EB_example.m```) for calling the main function is provided, along with synthetic example data (stored in the ```./example_data/ directory```). The same script can be used to run the function on any given network, simply by substituting the two input files:
1) a network in edge list form (space-separated, no header) 
	* columns: edge_start_node edge_end_node weight
	* weight represents the geographic length of the edge
2) a comma-separated file of 4 columns and N lines, N being the number of nodes in the network, and the columns being
	1) the x coordinate of the node (in a non-degree Coordinate Reference System)
	2) the y coordinate of the node (in a non-degree Coordinate Reference System)
	3) the population (i.e. number of pedestrians) assigned to the node
	4) the number of POIs assigned to the node

Note that the empirical sidewalk networks and node metadata used in [1] are stored in a persistent repository (https://www.doi.org/10.17605/OSF.IO/94TDC), and are correctly formatted and ready to use with a script like get_directed_EB_example.m

The script can be run from the MATLAB command line within the MATLAB GUI, and should execute very quickly (< 1 minute).

```
get_directed_EB_example
```

The progress of the algorithm will be reported in the command line window, starting from "Source=1" and ending at "Source=N" where N is the number of network nodes.

Otherwise, to call from a bash command line, run:

```
nohup matlab < get_directed_EB_example.m > output.txt &
```

where the progress of the algorithm will be written to ```./output.txt```. 

#### Output of the example script
The example script writes three files to a new folder ```./example_data/EB/```: ```eb_i.txt```, ```eb_j.txt```, and ```eb.txt```. These are formatted like the edge list input file: space-separated, three columns - edge_start_node edge_end_node edge_betweenness. ```eb_i.txt``` and ```eb_j.txt``` store the directed edge betweenness values of the upper and lower triangle of the output N x N array of directed edge betweenness values. ```eb.txt``` stores the undirected edge betweenness values of the network.

# References

[1] Rhoads, D., Solé-Ribalta, A., González, M. C., & Borge-Holthoefer, J. (2020). Planning for sustainable Open Streets in pandemic cities. arXiv preprint arXiv:2009.12548.


