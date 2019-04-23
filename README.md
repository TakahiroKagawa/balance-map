# balance-map
Matlab code for balance map analysis which evaluate walking balance.
Simulation results presented at ICRA2019 were produced by the matlab functions.
The simulations were built in matlab 2014b.

Matlab functions
BalanceMap.m plots boundary lines of fall regions. 
BalanceMapLimit.m calculates boundary lines of forward fall regions. So far, it will take a couple of minites to complete. There must be better algorithm on it.
TwoLinkModel_sim.m performs computer simulations of a swing movement with and without disturbance.
TwoLinkModel_comp.m performs computer simulations of a linear and a nonlinear compass gait model.
DefInitChi.m produces an initial state of the compass gait model for computer simulations.
StabilityBoundaryNegaitive_Model/StabilityBoundaryPositive_Model.m search the boundary state of fall regions.
TwoLinkFallDetection.m 
