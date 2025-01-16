%-------------------------------------------------------------------------------
%  [function]  Residual function: Finite elements (tents) + Galerking weights
%-------------------------------------------------------------------------------
function err= residuals_FE( theta,n, xg,wg,nq, kg, zg,nz,Pi, ...
                            beta,sigma,eta,alpha,delta )

    % unpack coefficients
    theta=  reshape(theta,n,nz);
    
    % vectorize for convenience
    k=  xg(:);                                      % capital sub-interval nodes
    w=  wg(:);                                      % sub-interval weights
    nq= nq*(n-1);                                   % number nodes

    % consumption {t}
    psi_q=  tent_basis( k,nq, kg,n );               % basis over quadrature nodes
    Chat=   psi_q* theta;                           % policy function
    % get capital policy function
    [~,Kp]= policy_functions( k,zg, Chat, sigma,eta,alpha,delta );

    
    % compute residuals
    wres= ones(n,nz);
    % loop over productivity
    for iz= 1:nz
        % current variables: over (k;z)
        c=      Chat(:,iz);                         % c{t}
        kp=     Kp(:,iz);                           % k{t+1}
        
        % future variables: over (k'(k;z),z')
        % consumption {t+1}
        Cp=     tent_basis(kp,nq, kg,n)*theta;
        % return on capital {t+1}
        [~,~,Rp]= policy_functions( kp,zg, Cp, sigma,eta,alpha,delta );
    
        % EE residuals
        res=    1  -  beta* c.^(sigma) .*(Rp .*Cp.^(-sigma)) *Pi(iz,:)';
        % weighted residuals
        res_i=      psi_q .*res;                    % Galerkin weighting
        wres(:,iz)= sum( res_i .* w )';             % integrate over capital
    end

    % vectorize residuals
    err= wres(:);
end
