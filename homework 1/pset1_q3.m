clear variables; close all; clc;

%===============================================================================
%  Q3. Optimization
%===============================================================================
% for printing figures
figdir=     'figures';
print_fig=  1;


% parameters
a=  100;
b=  1;

% initial guess
x0a= [0.25, 3]';
x0b= [1.75, -0.15]';


% routine parameters
max_iters= 1e5;
tol_x=  1e-8;
tol_g=  1e-5;

% gradient & Hessian functions
gradient=   @(x) gx(x,a,b);
Hessian=    @(x) Hx(x,a);

% benchmark: Matlab canned routine
tic;
f= @(x) a*( x(2)-x(1)^2 )^2 + (b-x(1))^2;
[xs0,f0]= fminsearch(f,x0a);
toc;

% I. Newton-Rhapson
tic;
[xs1,iter1,xk1,g1]= newton_rhapson( x0a, gradient,Hessian, max_iters,tol_x,tol_g );
toc;

% II. BFGS
tic;
H0= eye(length(x0a));                                    % initialize gradient
[xs2,iter2,g2]= BFGS( x0a,H0, gradient, max_iters,tol_x,tol_g );
toc;

% III. Steepest descent
alpha= 1e-2;
tic;
[xs3,iter3,g3]= steepest_descent( x0a, gradient, alpha, max_iters,tol_x,tol_g );
toc;

% IV. Conjugate descent
alpha= 1e-4;
tic;
[xs4,iter4,g4]= conjugate_descent( x0a, gradient, alpha, max_iters,tol_x,tol_g );
toc;


%-------------------------------------------------------------------------------
% figures
%-------------------------------------------------------------------------------
% graoh specs: (font,x,y)= (pt,in,in) 
fontsize= 9;  xsize= 4;  ysize= 3;
graph_specs(xsize,ysize,fontsize);

% objective
n= 1000;
xg= linspace(-1,2,n);
yg= linspace(-1,3,n);
F= zeros(n,n);
for i=1:n, for j=1:n, F(i,j)= f([xg(i) yg(j)]); end; end

% plot function
figure(1)
tiledlayout('flow');
nexttile;
hold on;
mesh(xg,yg,F');
view(-50,35);
p1= scatter3(xs0(1),xs0(2),f0, 'r', 'filled');
p2= scatter3(x0a(1),x0a(2),f(x0a), 'k', 'filled');
scatter3(x0b(1),x0b(2),f(x0b), 'k', 'filled');
hold off;
xlabel('x'); ylabel('y');
legend([p2,p1],'guesses','minimum', 'fontsize',fontsize );
lgd.Layout.Tile= 'south';
% [print]
if (print_fig==1), print(sprintf('%s\\fig_q3_function',figdir),'-depsc2'); end





%===============================================================================
%  Functions
%===============================================================================

%-------------------------------------------------------------------------------
%  gradient
%-------------------------------------------------------------------------------
function [df]= gx(xx,a,b)
    % unpack variables
    x= xx(1);  y= xx(2);
    % gradient
    df= zeros(2,1);
    df(1)=  4*a*x^3 + (2-4*a*y)*x - 2*b;
    df(2)=  2*a*(y-x^2);
end

%-------------------------------------------------------------------------------
%  Hessian
%-------------------------------------------------------------------------------
function [H]= Hx(xx,a)
    % unpack variables
    x= xx(1);  y= xx(2);
    % hessian
    H= zeros(2,2);
    H(1,1)= 12*a*x^2 + 2 - 4*a*y;
    H(2,2)= 2*a;
    H(1,2)= -4*a*x;
    H(2,1)= H(1,2);
end

%-------------------------------------------------------------------------------
%  Newton-Raphson method
%-------------------------------------------------------------------------------
function [x1,iter,xk,g0,H0]= newton_rhapson( x0, gradient,Hessian, max_iters,tol_x,tol_g )

    % preallocate path
    xk= zeros(max_iters,2);
    xk(1,:)= x0;

    % run Newton-Raphson routine
    for iter= 1:max_iters

        % gradient & Hessian
        g0= gradient(x0);
        H0= Hessian(x0);

        % update x
        x1= x0 - H0\g0;
        xk(iter+1,:)= x1;

        % error
        err_x= max(abs(x0-x1));
        err_g= norm(g0);
        x0= x1;
        % convergence criteria
        if ( err_x<tol_x || err_g<tol_g ); break; end
    end

    % store path
    xk= xk(1:iter+1,:);
end

%-------------------------------------------------------------------------------
%  Broyden–Fletcher–Goldfarb–Shanno method
%-------------------------------------------------------------------------------
function [x1,iter,g0,H0]= BFGS( x0, H0, gradient, max_iters,tol_x,tol_g )

    % run BFGS routine
    for iter= 1:max_iters

        % gradient
        g0= gradient(x0);

        % update x
        x1= x0 - H0\g0;
    
        % update hessian
        g1= gradient(x1);                                       % new gradient
        s=  x1 - x0;                                            % step size
        y=  g1 - g0;                                            % gradient change
        H0= H0 + (y*y')/(y'*s) - (H0*(s*s')*H0')/(s'*H0*s);     % next Hessian
    
        % error
        err_x= max(abs(x0-x1));
        err_g= norm(g0);
        x0= x1;
        % convergence criteria
        if ( err_x<tol_x || err_g<tol_g ); break; end
    end
end

%-------------------------------------------------------------------------------
%  Steepest Descent
%-------------------------------------------------------------------------------
function [x1,iter,g0]= steepest_descent( x0, gradient, alpha, max_iters,tol_x,tol_g )

    % run steepest-descent routine
    j= 2;
    for iter= 1:max_iters

        % alpha adjustment: no line search, try some decay
        if (mod(iter,10^j)==0), alpha= max(alpha/10,1e-7); j=j+1; end

        % gradient
        g0= gradient(x0);
        % descent direction
        d=  -g0/norm(g0);

        % update x
        x1= x0 + alpha*d;
    
        % error
        err_x= max(abs(x0-x1));
        err_g= norm(g0);
        x0= x1;
        % convergence criteria
        if ( err_x<tol_x || err_g<tol_g ); break; end
    end
end

%-------------------------------------------------------------------------------
%  Conjugate Descent
%-------------------------------------------------------------------------------
function [x1,iter,g0]= conjugate_descent( x0, gradient, alpha, max_iters,tol_x,tol_g )

    % initialize
    g0= gradient(x0);
    d=  -g0;
    x1= x0 + alpha*d;

    % run conjugate-descent method
    for iter= 1:max_iters
    
        % compute new direction
        g1= gradient(x1);                                       % new gradient
        bb= (g1'*g1)/(g0'*g0);                                  % dampening factor
        d=  -g1 + bb*d;                                         % new direction
        g0= g1;                                                 % last gradient
        
        % update x
        x1= x0 + alpha*d;

        % error
        err_x= max(abs(x0-x1));
        err_g= norm(g0);
        x0= x1;
        % convergence criteria
        if ( err_x<tol_x || err_g<tol_g ); break; end
    end
end
