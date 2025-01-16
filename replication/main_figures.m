%===============================================================================
%  program: replication exercises
%===============================================================================
%   purpose:    replicates parts of Achdou, Han, Lasry, Lions, & Moll (restud2021),
%               particurlarly the numerical results for the stationary setup.
% 
%   [notes]:
%     - parameters aren't the same across figures, bounds are set
%       differently in each set
%===============================================================================
clear variables; close all; clc;

% for printing figures
fig_dir=            'figures';
print_fig=          0;
fontsize=           10;
xsize=              8;
ysize=              3;


% common parameters
rho=                0.05;

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




%--------------------------------------------------
%  Figure 1
%--------------------------------------------------
% parameters
% asset grid
amin=   -0.15;
amax=   4;
na=     2000;
% productivity
ny=     2;
yg=     [0.1 0.2];
lambda= [-1.5 1.5; 1.6 -1.6];
% iterest rate & preferences
r=      0.025;
sigma=  2;

% I. set up
% build asset grid
da=     (amax-amin)/(na-1);
ag=     amin + ((1:na)'-1)*da;
% build stochastic transition matrix
Lambda= build_stoch_transition( lambda,na );

% II. solve HJB
[~,cpf,spf]= HJB_poisson( ag,da,na, yg,ny,Lambda, ...
                          r, rho,sigma,  Delta,tol_vf,max_iters_vf );



% [figure]: policy functions
% set plot bounds
a0=     -0.2;
a1=     0.4;
% legend
lgd=    {'$c(a,y_1)$','$c(a,y_2)$', ...
         'orientation','vertical', 'fontsize',fontsize-1};
% do figure
graph_specs(xsize,ysize,fontsize)
figure(1);
tiledlayout('flow');
% consumption policy
nexttile;
plot(ag,cpf);
xlim([a0 a1]);
xline( amin,'--k','$\underline{a}$');
xlabel('wealth, $a$');
ylabel('$c(a,y)$');
title('Consumption Policy Function');
legend(lgd{:},'location','SE');
% savings policy
nexttile;
plot(ag,spf);
xlim([a0 a1]);
xline( amin,'--k','$\underline{a}$');
yline(0,'k--');
xlabel('wealth, $a$');
ylabel('$\dot{a}(a,y)$');
title('Saving Policy Function');
legend(lgd{:},'location','NE');

% [print]
if (print_fig==1)
    print(sprintf('%s\\figure1',fig_dir),'-depsc2');
    exportgraphics(gcf, sprintf('%s\\figure1.pdf',fig_dir), 'ContentType', 'vector');
end


%--------------------------------------------------
%  Figure 4
%--------------------------------------------------
% parameters
% asset grid
na=     500;
amin=   -0.02;
amax=   3;
% productivity
ny=     2;
yg=     [0.1 0.2];
lambda= [-0.5 0.5; 0.3 -0.3];
% interest rate & preferences
sigma=  1.2;
r=      0.035;
% useful constants
eta=    (rho-r)/sigma;

% I. set up
% build asset grid
da=     (amax-amin)/(na-1);
ag=     amin + ((1:na)'-1)*da;
% build stochastic transition matrix
Lambda= build_stoch_transition( lambda,na );

% II. solve HJB
[~,cpf,~, A_ay,iiF,ii0]= HJB_poisson( ag,da,na, yg,ny,Lambda, ...
                                      r, rho,sigma,  Delta,tol_vf,max_iters_vf );

% III. compute MPCs
do_correction=  1;                                  % interpolate MPC errors
t1=             2;                                  % time-interval size
nt=             41;                                 % time discretization
% compute MPCs
[dcpf,cMPC]= MPC( ag,cpf, da,na,ny, t1,nt,A_ay, do_correction,iiF,ii0 );


% [figure]: policy functions
% set plot bounds
a0=     -0.1;
a1=     1.5;
y0=     0;
y1=     1.25;
% legend
lgd=    {'MPC$(a,y_1)$','MPC$(a,y_2)$', ...
         'orientation','vertical', 'location','NE', 'fontsize',fontsize-1};
% do figure
graph_specs(xsize,ysize,fontsize);
figure(4);
tiledlayout('flow');
% cummulative MPC
nexttile;
hold on;
plot(ag(1:na-1),cMPC);
plot(ag(1:na-1), ones(na-1,1)*(eta+r)*t1, '--g' );
hold off;
xline( amin,'--k','$\underline{a}$');
xlim([a0 a1]);
xlabel('wealth, $a$');
ylabel('MPC$_{j\tau}(a)$');
title('Cumulative MPC: over $[0,\tau]$');
legend(lgd{:});
% instantaneous MPC
nexttile;
hold on;
plot(ag(1:na-1),dcpf);
plot(ag(1:na-1), ones(na-1,1)*(eta+r), '--g' );
hold off;
xline( amin,'--k','$\underline{a}$');
xlim([a0 a1]);
ylim([y0 y1]);
xlabel('wealth, $a$');
ylabel('$\partial_a c(a,y_j)$');
title('Instantaneous MPC');
legend(lgd{:});

% [print]
if (print_fig==1)
    print(sprintf('%s\\figure4',fig_dir),'-depsc2');
    exportgraphics(gcf, sprintf('%s\\figure4.pdf',fig_dir), 'ContentType', 'vector');
end


%--------------------------------------------------
%  Figure 6
%--------------------------------------------------
% parameters
% asset grid
amin=   -0.15;
amax=   4;
na=     2000;
% productivity
ny=     2;
yg=     [0.1 0.2];
lambda= [-1.5 1.5; 1.6 -1.6];
% interest rate & preferences
sigma=  2;
r=      0.025;

% I. set up
% build asset grid
da=     (amax-amin)/(na-1);
ag=     amin + ((1:na)'-1)*da;
% build stochastic transition matrix
Lambda= build_stoch_transition( lambda,na );

% II. get HJB soln
[~,~,spf, A_ay]= HJB_poisson( ag,da,na, yg,ny,Lambda, ...
                              r, rho,sigma,  Delta,tol_vf,max_iters_vf );
% III. get KFE soln
g= KFE_poisson( A_ay, da,na,ny, Delta,tol_dist,max_iters_dist );


% [figure]: savings & stationary distribution
% set plot bounds
a0=     -0.2;
a1=     0.4;
a0d=    a0;
a1d=    0.4;
y0d=    0;
y1d=    4.5;
% legend
lgd=    {'orientation','vertical','location','NE', 'fontsize',fontsize-1};
graph_specs(xsize,ysize,fontsize)
figure(6);
tiledlayout('flow');
% savings policy
nexttile;
plot(ag,spf);
xlim([a0 a1]);
xline(amin,'--k','$\underline{a}$');
yline(0,'k--');
xlabel('wealth, $a$');
ylabel('$\dot{a}(a,y)$');
title('Saving Policy Function');
legend('$s(a,y_1)$','$s(a,y_2)$',lgd{:});
% distribution
nexttile;
plot(ag,g);
xline(amin,'--k','$\underline{a}$');
xlim([a0d a1d]);
ylim([y0d y1d])
xlabel('assets');
ylabel('$g(a,y_j)$');
title('Stationary Densities');
legend('$g(a,y_1)$','$g(a,y_2)$',lgd{:});

% [print]
if (print_fig==1)
    print(sprintf('%s\\figure6',fig_dir),'-depsc2');
    exportgraphics(gcf, sprintf('%s\\figure6.pdf',fig_dir), 'ContentType', 'vector');
end




%--------------------------------------------------
%  Figure 7
%--------------------------------------------------
% parameters
r1=     0.04;
% use other parameters from Figure 6

% II. get HJB soln
[~,~,spf1, A_ay1]= HJB_poisson( ag,da,na, yg,ny,Lambda, ...
                                r1, rho,sigma,  Delta,tol_vf,max_iters_vf );
% III. get KFE soln
g1= KFE_poisson( A_ay1, da,na,ny, Delta,tol_dist,max_iters_dist );


% [figure]: savings & stationary distribution
% set plot bounds
a0=     -0.2;
a1=     0.4;
a0d=    a0;
a1d=    0.4;
y0d=    0;
y1d=    4.5;
% legend
lgd=    {'orientation','vertical','location','NE', 'fontsize',fontsize-1};
graph_specs(xsize,ysize,fontsize)
figure(7);
tl= tiledlayout('flow');
% savings policy
nexttile;
hold on;
plot(ag,spf);
set(gca,'ColorOrderIndex',1);
plot(ag,spf1, '--');
hold off;
xlim([a0 a1]);
xline(amin,'--k','$\underline{a}$');
yline(0,'k--');
xlabel('wealth, $a$');
ylabel('$\dot{a}(a,y)$');
title('Saving Policy Function');
legend('$s(a,y_1)$','$s(a,y_2)$',lgd{:});

% distribution
nexttile;
hold on;
plot(ag,g);
set(gca,'ColorOrderIndex',1);
plot(ag,g1, '--');
hold off;
xline(amin,'--k','$\underline{a}$');
xlim([a0d a1d]);
ylim([y0d y1d])
xlabel('assets');
ylabel('$g(a,y_j)$');
title('Stationary Densities');
legend('$g(a,y_1)$','$g(a,y_2)$',lgd{:});

% [print]
if (print_fig==1)
    print(sprintf('%s\\figure7_Huggett',fig_dir),'-depsc2');
    exportgraphics(gcf, sprintf('%s\\figure7_Huggett.pdf',fig_dir), 'ContentType', 'vector');
end





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