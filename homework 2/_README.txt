Homework 2 codes:


  [programs]
	- main:			executes the whole program, reproducing all the results.
				(solves questions 1, 2, and 3, and produces graphs for question 5).

  [functions]
    * relevant for all questions:
	- discretize_productivity.m:	discretizes a demeaned AR(1) process with Tauchen's (1986) method.
	- policy_functions.m:		computes the labor, future capital, and associated capital return as described
					in equations (7)-(9).

    * relevant for question 1:
	- guess_chebyshev.m:		produces a coefficient guess for the Chebyshev collocation problem. 		[Question 1]
	- chebyshev_poly.m:		evaluates the Chebyshev polynomials.						[Question 1]
	- residuals_chebyshev.m:	computes the collocation residuals, as described in equations (5) and (10).	[Question 1]	

    * relevant for question 2:
	- guess_FE.m:			produces a coefficient guess for the finite elements problem. 			[Question 2]
	- tent_basis.m:			evaluates the tent basis, as described in equation (13). 			[Question 2]
	- residuals_FE.m:		computes the collocation residuals, as described in equations (14) and (15).	[Question 2]
	- gauss_legendre_quadrature.m:	Gauss-Legendre numerical integration with 10 nodes.				[Question 2]

    * relevant for question 3:
	- perturbation.mod:		contains the Dynare model. 							[Question 3]
	- perturbed_pf.m:		builds the policy function from the Dynare output. 				[Question 3]

    * relevant for question 5:
	- euler_eqn_residuals.m		computes the Euler equation residuals, as described in equation (16). 		[Question 5]
	- simulate_bilinear.m		simultes the economy's states using bilinear interpolation on k'(k,z).		[Question 5]
	- gaph_specs.m			sets global plotting specifications. 						[Question 5]
