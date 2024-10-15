# ddsip
Dual Decomposition in Stochastic Mixed-Integer Linear Programming

C implementation using IBM ILOG CPLEX and ConicBundle by C. Helmberg

For a description see the included manual ddsip-man.pdf (above)

seems that the use of CPLEX infeasilibility proof (Farkas ray) is not determininistic.
So be aware that each run could give a different performance due to added cuts.
