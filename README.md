# NodeGen
This is NodeGen. It is my attempt at generating a map for Strategy Games by first generating a graph and then to place the player spawns on graph's edges. The Graph's edges would represent how players would move between two vertices while vertices represent a choice between moving left or right (in the case of a vertex with order 3) or in other cases where the vertex has a higher order, a choice between many directions. 

The rest of the map would then be generated by generating terrain inside the regions of the graph and thus chokepoints could be created between the regions surrounding an edge.

Unfortunately I didn't get very far with the experiment, but it was fun and educational nonetheless.

# Notes
The library contains
- Point.java
- Delaunay.java
- Angle.java
- CoreGeom.java
All of which can be found in my ComputationalGeometry repo found at https://github.com/HeinrichWizardKreuser/ComputationalGeometry

