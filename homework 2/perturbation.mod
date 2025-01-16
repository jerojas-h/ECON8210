//------------------------------------------------------------------------------
//  Standard RBC model + endogenous labor supply
//------------------------------------------------------------------------------

// load parameter values
load parameters.mat;
set_param_value('beta'  ,parameters(1));             // discount factor
set_param_value('sigma' ,parameters(2));             // CRRA consumption parameter
set_param_value('eta'   ,parameters(3));             // CRRA labor parameter
set_param_value('alpha' ,parameters(4));             // capital intensity
set_param_value('delta' ,parameters(5));             // depreciation
// shock parameters
set_param_value('phi'   ,parameters(6));             // persistence
set_param_value('sd_e'  ,parameters(7));             // volatility


//------------------------------------------------------------------------------
//  parameters
//------------------------------------------------------------------------------
parameters beta, sigma, eta, alpha, delta, phi, sd_e
           rho y_k k_l;

// useful constants
rho=    1/beta-1;
y_k=    (rho+delta)/alpha;
k_l=    y_k^( 1/(alpha-1) );


//------------------------------------------------------------------------------
//  variable declaration block
//------------------------------------------------------------------------------
var     c, l, k, y, z;
varexo  e;

//------------------------------------------------------------------------------
//  optimality conditions
//------------------------------------------------------------------------------
model;

// euler equation
c^(-sigma)=  beta* ( alpha*exp(z(+1))*(k(+1)/l(+1))^(alpha-1) + (1-delta) ) *c(+1)^(-sigma);

// labor intratemporal FOC
l^eta=  (1-alpha)* k(-1)^alpha *l^(-alpha) *c^(-sigma);

// budget constraint
c + k=  y + (1-delta)*k(-1);

// output
y=  exp(z) *k(-1)^alpha *l^(1-alpha);

// productivity
z=  phi*z(-1) + e;

end;

//------------------------------------------------------------------------------
//  shocks block
//------------------------------------------------------------------------------
shocks;
var e; stderr sd_e;
end;

// to run simulations
stoch_simul( order= 3, irf=0);


//------------------------------------------------------------------------------
//  deterministic steady states
//------------------------------------------------------------------------------
initval;

// steady states
c=    ( k_l^alpha - delta*k_l ) *( (1-alpha)*(k_l)^alpha )^(1/eta);
c=    c^( 1/( 1 + sigma/eta ));
l=    ( (1-alpha)* k_l^alpha * c^(-sigma) )^(1/eta);
k=    k_l*l;
y=    y_k*k;
e=    0;

end;

