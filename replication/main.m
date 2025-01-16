%===============================================================================
%  program: general equilibrium
%===============================================================================
%   purpose:    for computing general equilibrium and testing
%               parametrization
%===============================================================================
clear variables; close all; clc;

% parameters
rho=    0.05;
sigma=  2;
% asset grid
na=     3000;
amin=   -0.15;
amax=   15;
% value function iteration
Delta=              1000;
tol_vf=             1e-6;
max_iters_vf=       100;
% distribution iteration
tol_dist=           1e-8;
max_iters_dist=     100;
% bisection
tol_A=              1e-4;
max_iters_bisect=   30;
rlb=                0.01;
rub=                rho;

% exogenous income parameters
y=      0.15;
yfac=   0.05;
ny=     3;


%--------------------------------------------------
%  set up grids & stochastic transition matrix
%--------------------------------------------------
% build productivity grid
yg=     linspace( y-yfac, y+yfac, ny );
lambda= 1.55*ones(1,ny);

% stochastic transitions
% iniatialize
Lambda= ones(ny,1)*lambda;
for iy= 1:ny
    % replace diagonal elements
    Lambda(iy,iy)=  -sum(lambda( (1:ny)~=iy ));
end
Lambda= build_stoch_transition( Lambda,na );

% build asset grid
da=     (amax-amin)/(na-1);
ag=     amin + ((1:na)'-1)*da;


%--------------------------------------------------
%  solve general equilibrium
%--------------------------------------------------
fprintf('Solving equilibrium interest rate...\n');
r= bisect_irate( rlb,rub, ag,da,na, yg,ny,Lambda, ...
                 tol_A,max_iters_bisect, tol_vf,max_iters_vf, tol_dist,max_iters_dist, ...
                 Delta, rho,sigma );

%--------------------------------------------------
%  solve HJB & KFE
%--------------------------------------------------
% solve HJB
[vf,cpf,spf, A_ay,iiF,ii0,iiB]= HJB_poisson( ag,da,na, yg,ny,Lambda, ...
                                             r, rho,sigma,  Delta,tol_vf,max_iters_vf );
% solve KFE soln
[g,G]= KFE_poisson( A_ay, da,na,ny, Delta,tol_dist,max_iters_dist );





%-------------------------------------------------------------------------------
%  Figures
%-------------------------------------------------------------------------------
% get bounds for ploting
qa= 0.999;
amax_dist=  ag( sum(sum(G,2)<qa) );                 % upper bound for plots


% [figure]: agent's problem
graph_specs(8,6,12)
figure(1);
tl= tiledlayout('flow');
% consumption policy
nexttile;
plot(ag,cpf);
xlim([amin amax]);
xlabel('assets');
ylabel('consumption rate');
title('$c(a,y)$');
% savings policy
nexttile;
plot(ag,spf);
xlim([amin amax]);
yline(0);
xlabel('assets');
xlabel('saving rate');
title('$\dot{a}(a,y)$');
% value function
nexttile;
plot(ag,vf);
xlim([amin amax]);
xlabel('assets');
xlabel('value');
title('$v(a,y)$');
% distribution
nexttile;
plot(ag,G);
xlim([amin amax_dist]);
xlabel('assets');
ylabel('CDF');
title('$G(a,y)$');

figure(2);
plot(ag,g);
xlim([amin amax_dist]);
% ylim([0 .015]);
xlabel('assets');
ylabel('pdf');
title('$g(a,y)$');




%===============================================================================
%  Functions
%===============================================================================

%-------------------------------------------------------------------------------
%  [function]:  stochastic transition matrix
%-------------------------------------------------------------------------------
function [Lambda]= build_stoch_transition( Lambda, na )

    % map to (ay x ay) space
    Lambda= kron( Lambda, eye(na) );
    Lambda= sparse(Lambda);

end