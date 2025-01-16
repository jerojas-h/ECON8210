This repository contains the codes for problem set 1 of ECON8210.
	Date: 10/31/2024

All questions are self-contained in their own m-files:
	- pset1_q2.m:		solves question 2 - integration
	- pset1_q3.m:		solves question 3 - optimization
	- pset1_q4_q5.m:	solves question 4 & 5 - Pareto efficient allocations
						      & Competitive equilibria.
	- pset1_q6.m:		solves question 6, currently parts 2 and 3 only.
	- graph_specs.m:	sets global template for figures.
		

Additional documentation:

  - pset1_q2.m:		contains nested functions for the Newton-Coates routines and
			  Monte Carlo integration (not optimized for speed).

  - pset1_q3.m:		contains nested functions to evaluate the analytic gradient
			  and Hessian.
			contains nested functions for each method: (i) Newton-Rhapson,
			  (ii) BFGS, (iii) Steepest Descent, and (iv) Conjugate Descent.

  - pset1_q4_q5.m:	contains nested functions computing (i) the planner's resource
			  constraint errors and (ii) the decentralized equilibrium 
			  condition errors, both described in the pdf.
			questions 4 and 5 are solved using 'fsolve'.

  - pset_q6.m:		currently incomplete, described in the .pdf document

