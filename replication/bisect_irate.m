%-------------------------------------------------------------------------------
% [Function]:  interest rate bisection
%-------------------------------------------------------------------------------
%   Purpose:  computes the equilibrium interest rate.
% 
%   [notes]:
%     - currently runs for poisson processes
%-------------------------------------------------------------------------------
function [r,vf,cpf,spf, g,G, A_ay]= bisect_irate( rlb,rub, ag,da,na, yg,ny,Lambda, ...
    tol_A,max_iters_bisect, tol_vf,max_iters_vf, tol_dist,max_iters_dist, ...
    Delta, rho,sigma )

    % initialize bisection
    r=  (rlb+rub)/2;
    
    for iter_bisection= 1:max_iters_bisect
    
        % solve HJB & transition matrix
        [vf,cpf,spf, A_ay]= HJB_poisson( ag,da,na, yg,ny,Lambda, r, ...
                                         rho,sigma, Delta,tol_vf,max_iters_vf );
        %  solve KFE
        [g,G]= KFE_poisson( A_ay, da,na,ny, Delta,  tol_dist,max_iters_dist);
    
        %------------------------------------------
        %  compute excess demand
        %------------------------------------------
        % asset demand
        Da= kron(ones(ny,1),ag)' *g(:)*da;
        % excess demand
        err_A=  Da;
    
        % [print] progress
        fprintf('iter: %3d  |  r= %.2f%%  |  Z(a)= %.2e  \n', ...
            iter_bisection, 100*r, err_A );
        % check convergence
        if (abs(err_A)<tol_A), break; end
    
        % do bisection
        if (err_A<0)
            % increase interest rate
            r0=  r;
            r=   (r+rub)/2;   
            rlb= r0; 
        else
            % decrease interest rate
            r0=  r;
            r=   (r+rlb)/2;   
            rub= r0; 
        end
    
    end

end