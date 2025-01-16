clear variables; close all; clc;

%===============================================================================
%  Questions 4 & 5: set up
%===============================================================================
opt= optimset('Display','off');

% parameters
n=  30;
m=  30;

% preference parameters
alph=   linspace(1,m,m) + ones(1,m)*0;
alph=   alph/sum(alph);
omg=    linspace(-4,-2,m) +  ones(1,m)*(-2)*0;
% planner weights
lambdas= linspace(1,n,n).^2/n;
lambdas= lambdas/sum(lambdas);

% generate preference parameter
Alphas= zeros(n,m);
Omegas= zeros(n,m);
for i= 1:n
    for j= 1:m
        k= i-1+j;
        k= k-m*(k>m);
        Alphas(i,k)= alph(j);
        Omegas(i,k)= omg(j);
    end
end
% generate endowments
Endowments= abs(Omegas)./sum(abs(Omegas),2) .*(1:m);



%-------------------------------------------------------------------------------
%  Q4. Computing Pareto efficient allocations
%-------------------------------------------------------------------------------
% aggregate endowment
Ek=  sum(Endowments);

% solve planner's problem
X= zeros(k,1);
for k= 1:m
    x1k= Ek(k)*lambdas(1);
    % compute allocation x{i,k}
    X(1,k)= fsolve( @(x) planner_rc( x,k,n, lambdas,Alphas,Omegas, Ek(k) ), x1k, opt );
    % compute remaining allocations
    for i= 2:n
        X(i,k)= ( lambdas(i)/lambdas(1)* Alphas(i,k)/Alphas(1,1) )^(-1/Omegas(i,k) ) ...
                * X(1,k)^(Omegas(1,1)/Omegas(i,k));
    end
end



%-------------------------------------------------------------------------------
%  Q5. Computing Pareto efficient allocations
%-------------------------------------------------------------------------------
% initialize
pk= ones(m-1,1)/(n*m);
mu= ones(n,1)/(n*m);
x0= [pk;mu];

% solve decentralized equilibrium
[xs,err]= fsolve( @(x) err_decentralized_equilibriun(x, n,m, Alphas,Omegas,Endowments), x0 );

% check Walras law
[~,Demand]= err_decentralized_equilibriun( xs,n,m, Alphas,Omegas,Endowments );
z1= sum(Demand(:,1)-Endowments(:,1));
fprintf('Market 1 excess demand: %.2e\n', z1 );

% fprintf('\nDemand:\n')
% disp(Demand);



%===============================================================================
%  Functions
%===============================================================================

%-------------------------------------------------------------------------------
%  planner's error
%-------------------------------------------------------------------------------
function [err]= planner_rc( x1k,k,n, lambdas,Alphas,Omegas, Ek )
    
    % summation term
    X= 0;
    for i= 1:n
        % add allocation x{i,k}
        X= X + ( lambdas(i)/lambdas(1)* Alphas(i,k)/Alphas(1,1) )^(-1/Omegas(i,k) ) ...
               * x1k^(Omegas(1,1)/Omegas(i,k));
    end
    % FOC error
    err= X - Ek; 

end

%-------------------------------------------------------------------------------
%  competitive equilibrium error
%-------------------------------------------------------------------------------
function [Err,Demand]= err_decentralized_equilibriun( x, n,m, ...
                                                      Alphas,Omegas,Endowments )
    % unpack variables
    price_k=    x(1:m-1);
    mu_i=       x(m-1+(1:n)); 
    % add normalized price
    price_k=    [ 1; price_k ];


    % demand
    Demand= zeros(n,m);
    for i= 1:n
        for k= 1:m
            Demand(i,k)= ( Alphas(i,k)* mu_i(i)/price_k(k) )^(-1/Omegas(i,k)); 
        end
    end

    % I. Excess demand
    Zk= sum( Demand-Endowments, 1)';
    Zk= Zk(2:m);

    % II. Budget constraint errors
    Zi= (Demand - Endowments)*price_k;

    % collect errors
    Err= [Zk;Zi];
end


