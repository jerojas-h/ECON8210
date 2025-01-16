Replication exercises:	this folder contains codes to reproduce select graphs from Achdou Han Lasry Lions Moll (restud2022). 


  [programs]
	- main:				for testing and solving the general equilibrium.
	- main_figures:			reproduces the selected figures.

  [functions]
	- HJB_poisson.m:		solves the HJB for an n-state Poisson income process using an implicit method
					with an upwind scheme.
	- KFE_poisson.m:		solves the KFE by iterating the distribution with the transition matrix 
					obtained from HJB_poisson.m.
	- MPC.m:			computes the instant and cumulative marginal propensities to consume.
	- bisect_irate.m:		solves the general equilibrium by performing an interest rate bisection.
	- gaph_specs.m			sets global plotting specifications.