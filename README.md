# MatrixFunctionRecycling
MATLAB code which demonstrates three differnt implementations of a new augmented Krylov subspace algorithm rFOM2 
allowing for krylov subspace recycling for matrix functions.

This code is available as part of the preprint "Krylov subspace recycling for matrix functions
- Liam Burke, Andreas Frommer, Gustavo Ramirez Hidalgo and Kirk M. Soodhalter (2022)".

The reader is encouraged to read the preprint before using the code. 

The code contains the following main run scripts

- quad_test.m     Tests the three implementations of rFOM2 as an augmented Krylov subspace method.
- recycle_test.m  Tests the three implementations of rFOM2 as a recycle method.

There is also two other (less important) run scripts

- U_test.m        Tests the quality of U as an eigenvector approximation.
- eigs_test.m     Looks at the minimum eigenvalue of each matrix in the sequence of problems. 

Any questions or comments on this code can be directed to Liam Burke (burkel8@tcd.ie). 
