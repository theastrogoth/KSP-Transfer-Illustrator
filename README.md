# KSP-Transfer-Illustrator
Work-in-progress KSP Transfer Illustrator, implemented in MATLAB

Use "testScript.m" to generate a porkchop plot and illustrate trajectories for the parameters set in the top input section.

Classes:
  Orbit - defines a Keplerian orbit via its 6 Keplerian orbital elements, as well as the body around which it rotates.
  Body - defines a celestial body by its name, size, gravity parameter, orbit, and color.
  Transfer - defines a solution to Lambert's problem for a single choice of start time and flight time, between two bodies.
  TransferTable - defines transfers across many start and end times in a table. 
  
  Functions:
    Transform3D - applies rotations on a vector to represent it by a new set of input bases.
    solveKepler - solves Kepler's equation via Newton's method for both elliptical and hyperbolic orbits.
    solveLambert - solves lambert's problem for a starting orbit around one body to an ending orbit around another, for a given start time and flight duration. 
