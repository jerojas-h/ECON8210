%-------------------------------------------------------------------------------
% [program] part 1-3 & 5
%-------------------------------------------------------------------------------
%   Purpose: solves questions 1, 2, 3, & 5.
% 
%   [notes]: 
%     - I pulled the quadrature nodes from numerical methods in C for n=10.
%-------------------------------------------------------------------------------
clear variables; close all; clc;

% for printing figures
fig_dir=    'figures';
print_fig=  0;
fontsize=   12;



% define parameters
beta=   0.97;
sigma=  1;
eta=    1;
alpha=  0.33;
delta=  0.1;

% state-space & basis
p1=     6;                                          % Chebyshev order
n2=     8;                                          % number of elements
kssfac= 0.25;
dnk=    1001;

% quadrature
nz=     3;
mz=     1;
sig_e=  0.007;
phi_z=  0.95;

% fsolve options
opt0= optimset('Display','off');
opt= optimset('Display','Iter','TolFun',10^(-15),'TolX',10^(-15));


%--------------------------------------
%  deterministic steady-states
%--------------------------------------
% useful constants
rho=    1/beta-1;
y_k=    (rho+delta)/alpha;
k_l=    y_k^( 1/(alpha-1) );
% steady-states
css=    ( k_l^alpha - delta*k_l ) *( (1-alpha)*(k_l)^alpha )^(1/eta);
css=    css^( 1/( 1 + sigma/eta ));
lss=    ( (1-alpha)* k_l^alpha * css^(-sigma) )^(1/eta);
kss=    k_l*lss;
yss=    y_k*kss;
xss=    [css lss kss yss ]';

%--------------------------------------
%  discretize state space
%--------------------------------------
% Tauchen's method
[zg,Pi]= discretize_productivity( nz, phi_z,sig_e, mz );

% dense capital grid
kmin=   (1-kssfac)*kss;
kmax=   (1+kssfac)*kss;
dkg=    linspace(kmin,kmax,dnk)';





%-------------------------------------------------------------------------------
%  Q1. Chebyshev collocation
%-------------------------------------------------------------------------------
fprintf('--------------------------------------------------------------------\n')
fprintf('Solving Chebyshev collocation...\n')
fprintf('--------------------------------------------------------------------')

% Chebyshev roots
Troots= -cos( (2*(1:p1)'-1)/(2*p1) *pi );             % order-p polynomial
% capital grid
kg_1=   kmin + (kmax-kmin)/2*( Troots+1 );          % capital grid

% guess coefficients
theta0= guess_chebyshev( kg_1,p1, zg,nz, beta,alpha, opt0 );
% do collocation
tic;
theta=  fsolve( @(th) residuals_chebyshev( th,p1, kg_1, zg,nz,Pi, ...
                                           beta,sigma,eta,alpha,delta ), theta0, opt );
toc;

% get coefficients
theta_1=        reshape(theta,p1,nz);
% compute policy functions
cpf_1=          chebyshev_poly(dkg,dnk,p1)*theta_1;
[lpf_1,kpf_1]=  policy_functions( dkg,zg, cpf_1, sigma,eta,alpha,delta );


%------------------------------------------
%  simulate ergodic distribution
%------------------------------------------
fprintf('\n\nSimulating economy...\n');
% parameters
T=  1*1e5;
T0= 1e4;
% set seed
rng(0115);

% get productivity innovations
T1=  T+T0;
et= sig_e*randn(T1,1);

% simulate states
tic;
[kt_bilinear,zt_bilinear]= simulate_bilinear( dkg,zg,kpf_1, kss,et,phi_z, T0,T1 );
toc;





%-------------------------------------------------------------------------------
%  Q2. Finite elements: tent-basis + Galerkin weighting
%-------------------------------------------------------------------------------
% number of elements
do_refined_guess= 1;

% build capital grid
kmin=   (1-kssfac)*kss;
kmax=   (1+kssfac)*kss;
kg_2=   kmin + (kmax-kmin)/(n2-1)* (0:n2-1)';

% set up quadrature for capital
nq= 10;                                             % #quadrature nodes
xg= zeros(nq,n2-1);
wg= zeros(nq,n2-1);
for i= 1:n2-1
    % get Gauss-Legendre nodes & weights for capital subintervals
    [xg(:,i),wg(:,i)]= gauss_legendre_quadrature(kg_2(i),kg_2(i+1));
end


% get inial guess: collocation + B1-splines
fprintf('--------------------------------------------------------------------\n')
fprintf('Guess: collocation w/ B1-splines...\n')
fprintf('--------------------------------------------------------------------')
tic;
theta= guess_FE( kg_2,n2, zg,nz,Pi, beta,sigma,eta,alpha,delta, do_refined_guess,opt );
th0= theta(:);
toc;

% solve FEM w/ Galerkin weights
fprintf('\n\n--------------------------------------------------------------------\n')
fprintf('FEM: tents + Galerkin weights...\n')
fprintf('--------------------------------------------------------------------')
tic;
theta= fsolve( @(theta) residuals_FE( theta,n2, xg,wg,nq, kg_2, zg,nz,Pi, ...
                                      beta,sigma,eta,alpha,delta ), th0, opt );
toc;
% store coefficients over (k,z)
theta_2=        reshape(theta,n2,nz);

% policy functions
cpf_2=          tent_basis( dkg,dnk, kg_2,n2 )*theta_2;
[lpf_2,kpf_2]=  policy_functions( dkg,zg, cpf_2, sigma,eta,alpha,delta );





%-------------------------------------------------------------------------------
%  Q3.  Perturbation
%-------------------------------------------------------------------------------
% store parameters to pass in Dynare
parameters= [beta sigma eta alpha delta phi_z sig_e];
save('parameters.mat','parameters');                    

% do perturbation with dynare
dynare perturbation.mod;

% check steady-states to make sure parametrization is the same
fprintf('\n\nCheck steady-states: \n');
disp([ xss oo_.steady_state(1:4) ]);

% get decision rules & approximation order
DR=             oo_.dr;
approx_order=   options_.order;

% state grids in steady-state deviations
% capital
dK= dkg-kss;
% productivity
dZ=     zg;

% switch the default ordering to the desired one for plotting
% (this is a lazy implementation)
pf_nam= {'c','l','k'};
pf_idx= zeros(length(pf_nam),1);
for i= 1:length(pf_nam)
    idx=        find(contains( oo_.var_list, pf_nam{i} ));
    pf_idx(i)=  find(oo_.dr.order_var==idx);
end

% compute policy functions
cpf_3= perturbed_pf( dK,dZ, DR,pf_idx(1), approx_order );
lpf_3= perturbed_pf( dK,dZ, DR,pf_idx(2), approx_order );
kpf_3= perturbed_pf( dK,dZ, DR,pf_idx(3), approx_order );






%===============================================================================
%  Figures
%===============================================================================
close all;
col= colororder;                                % get default colors

% legend
lgd=    {'Chebyshev collocation, $p=6$', 'Finite Elements + Galerkin, $n=8$', ...
         'Perturbation, order 3.', ...
         'fontsize',fontsize-1};


%------------------------------------------
%  Chebyshev vs FEM
%------------------------------------------
% [figure]
graph_specs(10,12,fontsize)
figure(1);
tiledlayout(3,2);
% consumption
nexttile;
hold on;
plot(dkg,cpf_1, 'Color',col(1,:) );
plot(dkg,cpf_2, '--', 'Color',col(2,:));
xline(kss,'--');
hold off;
xlim([kmin kmax]);
xlabel('capital');
ylabel('consumption');
title('$\hat{c}(k,z;\theta)$');
% labor
nexttile;
hold on;
plot(dkg,lpf_1, 'Color',col(1,:) );
plot(dkg,lpf_2, '--', 'Color',col(2,:));
xline(kss,'--');
hold off;
xlim([kmin kmax]);
xlabel('capital');
ylabel('labor');
title('$l(k,z;\theta)$');
% future capital
nexttile;
hold on;
plot(dkg,kpf_1, 'Color',col(1,:) );
g2= plot(dkg,kpf_2, '--', 'Color',col(2,:));
plot([kmin kmax],[kmin kmax], '--k');
xline(kss,'--');
hold off;
xlim([kmin kmax]);
ylim([ min(kpf_1(:)) , max(kpf_1(:)) ]);
xlabel('capital');
ylabel('future capital');
title('$k^\prime(k,z;\theta)$');


%------------------------------------------
%  Chebyshev vs perturbation
%------------------------------------------
% consumption
nexttile;
hold on;
plot(dkg,cpf_1, 'Color',col(1,:) );
plot(dkg,cpf_3, '--', 'Color',col(3,:));
xline(kss,'--');
hold off;
xlim([kmin kmax]);
xlabel('capital');
ylabel('consumption');
title('$\hat{c}(k,z;\theta)$');
% labor
nexttile;
hold on;
plot(dkg,lpf_1, 'Color',col(1,:) );
plot(dkg,lpf_3, '--', 'Color',col(3,:));
xline(kss,'--');
hold off;
xlim([kmin kmax]);
xlabel('capital');
ylabel('labor');
title('$l(k,z;\theta)$');
% future capital
nexttile;
hold on;
g1= plot(dkg,kpf_1, 'Color',col(1,:) );
g3= plot(dkg,kpf_3, '--', 'Color',col(3,:));
xline(kss,'--');
plot([kmin kmax],[kmin kmax], '--k');
hold off;
xlim([kmin kmax]);
ylim([ min(kpf_1(:)) , max(kpf_1(:)) ]);
xlabel('capital');
ylabel('future capital');
title('$k^\prime(k,z;\theta)$');

% legend
lgd_tl= legend( [g1(1) g2(1) g3(1)], lgd{:} );
lgd_tl.Layout.Tile= 'south';

% [print]
if (print_fig==1)
    print(sprintf('%s\\policy_functions',fig_dir),'-depsc2');
    exportgraphics(gcf, sprintf('%s\\policy_functions.pdf',fig_dir), 'ContentType', 'vector');
end



%------------------------------------------
%  euler equation residuals
%------------------------------------------
% I. Chebyshev collocation
method1.method= 1;
method1.theta=  theta_1;
method1.p=      p1;
% compute residuals
ee_res_1=   euler_eqn_residuals( cpf_1,kpf_1, method1, zg,nz, Pi, ...
                                 beta,sigma,eta,alpha,delta );
% II. finite elements method
method2.method= 2;
method2.theta=  theta_2;
method2.n=      n2;
method2.kg=     kg_2;
% compute residuals
ee_res_2=   euler_eqn_residuals( cpf_2,kpf_2, method2, zg,nz, Pi, ...
                                 beta,sigma,eta,alpha,delta );

% III. perturbation
method3.method= 3;
method3.kss=    kss;
method3.DR=     DR;
method3.idx_c=  pf_idx(1);
method3.order=  approx_order;
% compute residuals
ee_res_3=   euler_eqn_residuals( cpf_3,kpf_3, method3, zg,nz, Pi, ...
                                 beta,sigma,eta,alpha,delta );


% [figure]: Euler equation residuals
graph_specs(8,6,fontsize)
figure(2);
tiledlayout('flow');
nexttile;
% left axis
yyaxis left;
histogram( kt_bilinear, 'Normalization','probability', ...
           'FaceAlpha',0.25, 'EdgeAlpha', 0 );
xlim([kmin kmax]);
ylabel('density');
% right axis
yyaxis right;
hold on;
set(gca,'ColorOrderIndex',1);
g1= plot(dkg,log10(abs(ee_res_1)), 'Color',col(1,:), 'Marker', 'o');
g2= plot(dkg,log10(abs(ee_res_2)), 'Color',col(2,:), 'Marker', 'diamond');
g3= plot(dkg,log10(abs(ee_res_3)), 'Color',col(3,:), 'Marker', '^');
hold off;
xlim([kmin kmax]);
ylabel('$\mathcal{O}$(residuals)');
xlabel('capital');
% legend
lgd_tl= legend( [g1(1) g2(1) g3(1)], lgd{:} );
lgd_tl.Layout.Tile= 'south';

% right axis color
ax= gca;
ax.YAxis(2).Color= 'k';

% [print]
if (print_fig==1)
    print(sprintf('%s\\euler_equation_residuals',fig_dir),'-depsc2');
    exportgraphics(gcf, sprintf('%s\\euler_equation_residuals.pdf',fig_dir), 'ContentType', 'vector');
end