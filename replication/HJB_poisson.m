%-------------------------------------------------------------------------------
% [Function]:  Hamilton-Jacobi-Bellman equation w/ Poisson jumps
%-------------------------------------------------------------------------------
%   Purpose:  solves the stationary HBJ eqn for CRRA utility & endowments
%             following a Poisson process.
% 
%   [notes]:
%     - solved with upwind scheme: current scheme works if the value
%     function is concave.
%-------------------------------------------------------------------------------
function [vf,cpf,spf, A_ay, iiF,ii0,iiB]= HJB_poisson( ag,da,na, yg,ny,Lambda, ...
    r, rho,sigma,  Delta,tol,max_iters )

    % get parameters
    amax=   ag(na);
    amin=   ag(1);


    % initialize: value function guess
    v0=     (r*ag + yg).^(1-sigma)/(1-sigma)  /rho;     
    
    % iterate value function
    for iter= 1:max_iters
    
        %---------------------------------------
        % I. set up upwind-scheme directions
        %---------------------------------------
        % derivative approximation: finite differences
        dv=     ( v0(2:na,:) - v0(1:na-1,:) )/da;
        % forward difference
        dv_max= ( r*amax + yg ).^(-sigma);              % impose derivative on bound
        dvF=    [ dv; dv_max ];
        % backward difference
        dv_min= ( r*amin + yg ).^(-sigma);              % boundary condition: borrowing constraint
        dvB=    [ dv_min; dv ];
        
        % implied consumption possibilities
        cF=     dvF.^(-1/sigma);                        % FOC w/ forward difference
        cB=     dvB.^(-1/sigma);                        % FOC w/ backward difference
        c0=     r*ag + yg;                              % zero savings rate
        
        % savings functions
        sF=     r*ag + yg - cF;                         % w/ forward difference
        sB=     r*ag + yg - cB;                         % w/ backward difference
        % state regions: finite-difference direction
        iiF=    (sF>0);                                 % always positive savings
        iiB=    (sB<0);                                 % always negative savings
        ii0=    (1-iiF-iiB);                            % no savings (if concave value function)
        
        %---------------------------------------
        % II. update value function
        %---------------------------------------
        % instantaneous utility
        cpf=    iiF.*cF + ii0.*c0 + iiB.*cB;            % consumption policy
        u=      cpf.^(1-sigma) /(1-sigma);              % utility   
        
        % build transition matrix
        % 1. transition over a
        splus=  sF.*iiF /da;                            % forward-difference term
        splus=  splus(:);                               % vectorize
        sminus= sB.*iiB /da;                            % backward-difference term
        sminus= sminus(:);                              % vectorize
        % transition over a
        A_a=    spdiags( [-sminus (sminus-splus) splus], [1 0 -1], na*ny,na*ny )';
        % 2. transition over (a,y)
        A_ay=   A_a + Lambda;                           % add Poisson jumps
    
        
        % get value function update
        A_lhs=  ( 1/Delta + rho )*speye(na*ny) - A_ay;         
        u_vec=  u(:);                                   % vectorize u(c(a,y))
        v0_vec= v0(:);                                  % vectorize v(a,y)
        v1_vec= A_lhs\( u_vec + v0_vec/Delta );         % next value function
        
        % update value
        v0= reshape(v1_vec,na,ny);
        
        % check for convergence
        err_vf= max(abs(v0_vec-v1_vec));
        if (err_vf<tol), break; end
    
    end

    % [flags]
    if (iter==max_iters), fprintf('ERROR: max HJB iterations exceeded\n'); end

    % [output]
    vf=     v0;                                         % value function
    spf=    r*ag + yg - cpf;                            % savings policy function

end