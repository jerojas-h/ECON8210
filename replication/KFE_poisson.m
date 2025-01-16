%-------------------------------------------------------------------------------
% [Function]:  Kolmogrov Forward equation w/ Poisson jumps
%-------------------------------------------------------------------------------
%   Purpose:  solves the stationary KFE eqn for endowments following a Poisson
%             process.
% 
%   [notes]:
%     - the distribution is solved iteratively.
%-------------------------------------------------------------------------------
function [g,G]= KFE_poisson( A_ay, da,na,ny, Delta,  tol,max_iters)  

    % initialize
    Ap=     ( speye(na*ny) - Delta*A_ay' );         % set up transition matrix
    g0_vec= ones(na*ny,1)/(na*ny);                  % guess density
    

    % iterate forward in time
    for iter= 1:max_iters
        % compute next density
        g1_vec= Ap\g0_vec;

        % check convergence
        err_g=  max(abs(g0_vec-g1_vec));
        if (err_g<tol), break; end

        % update
        g0_vec= g1_vec;
    end
    

    % normalize density
    measure=    da* ones(1,na*ny)*g0_vec;           % implied measure
    g_vec=      g0_vec/measure;                     % for integrating to one
    
    % distributions
    g=  reshape( g_vec, na,ny);                     % pdf over (a,y)
    G=  cumsum(g*da,1);                             % CDF over (a,y)

end