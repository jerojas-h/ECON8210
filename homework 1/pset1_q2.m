clear variables; close all; clc;

%===============================================================================
%  Q2. Integration
%===============================================================================
% for printing figures
figdir=     'figures';
print_fig=  1;


% parameters
T=  100;
rho=    0.04;
lambda= 0.02;

% integrand
uc= @(c) -exp(-c);
fx= @(x) exp(-rho*x).*uc( 1-exp(-lambda*x) );


% grid sizes: newton coates
m_nc=   100;
ng_nc=  logspace(1,4,m_nc);
ng_nc=  floor(ng_nc/2)*2;
% grid sizes: monte carlo
m_mc=   100;
ng_mc=  logspace(1,8,m_mc);
ng_mc=  floor(ng_mc/2)*2;

% integral bounds
a= 0;  b= T;

% I. Newton-Coates integrals
IInc= zeros(m_nc,3);
tic;
for i= 1:m_nc
    nx= ng_nc(i);
    IInc(i,:)= integrate_newton_cotes(fx,a,b,nx);
end
toc;

% II. Monte Carlo integral
IImc= zeros(m_mc,1);
tic;
for i= 1:m_mc
    nx= ng_mc(i);
    IImc(i)= integrate_monte_carlo(fx,a,b,nx);
end
toc;


% comparison point
I=  integral(fx,a,b);


%-------------------------------------------------------------------------------
% figures
%-------------------------------------------------------------------------------
% graoh specs: (font,x,y)= (pt,in,in) 
fontsize= 10;  xsize= 4;  ysize= 3;
graph_specs(xsize,ysize,fontsize);

% errors
err_nc= log10(abs(IInc-I));
err_mc= log10(abs(IImc-I));

% newton-coates
figure(1);
tiledlayout('flow');
nexttile;
plot( ng_nc, err_nc );
xlabel('number of nodes');
ylabel('order of error' );
lgd= legend('Midpoint','Trapezoidal','Simpson', 'fontsize',fontsize );
lgd.Layout.Tile= 'south';
% [print]
if (print_fig==1), print(sprintf('%s\\fig_q2_newton_coates',figdir),'-depsc2'); end

% monte carlo
figure(2);
tiledlayout('flow');
nexttile;
plot( ng_mc, err_mc );
xlabel('number of nodes');
ylabel('order of error' );
% [print]
if (print_fig==1), print(sprintf('%s\\fig_q2_monte_carlo',figdir),'-depsc2'); end





%===============================================================================
%  Functions
%===============================================================================

%-------------------------------------------------------------------------------
%  Newton-Cotes integration rules
%-------------------------------------------------------------------------------
function [II]= integrate_newton_cotes(fx,a,b,n)

    % step size
    h=  (b-a)/n;
    
    % I. midpoint rule
    I1= h*sum(fx( a + (0.5:n)'*h ));
    
    % II. trapezoidal rule
    I2= h*sum(fx( a + (1:n-1)'*h )) + ...
        h/2*( fx(a) + fx(b) );
    
    % III. simpson's rule
    I3= h*4/3*sum(fx( a + (1:2:n-1)'*h )) + ...
        h*2/3*sum(fx( a + (2:2:n-1)'*h )) + ...
        h/3*( fx(a) + fx(b) );
    
    % [output]
    II= [I1 I2 I3]';
end

%-------------------------------------------------------------------------------
%  Monte Carlo integration
%-------------------------------------------------------------------------------
function [I]= integrate_monte_carlo(fx,a,b,n)
    % set seed for comparison
    rng(1031);
    
    % monte carlo integration
    m=  n;
    x=  a + (b-a)*rand(m,1);
    I=  (b-a)*sum(fx(x))/m;
end

